
const SCALING_FACTOR = 1_000_000

# use specialty type for "read counts"
# with basic function: 'sum' that can be
# overridden for more complex count bias models
#const ReadCount = JointBiasCounter #DefaultCounter
const DEFAULT_DEF   = (v = zero(DefaultCounter); v.isadjusted = true; v)
const DEFAULT_JOINT = (v = zero(JointBiasCounter); v.isadjusted = true; v)

const DEF_READCOUNT = 0.0
const DEF_READVALUE = 1.0

default(::Type{DefaultCounter}) = DEFAULT_DEF
default(::Type{JointBiasCounter}) = DEFAULT_JOINT

# Use this abstraction to store SGAlignPaths before
# quantification or assignment
abstract type SGAlignContainer end

struct SGAlignSingle <: SGAlignContainer
   fwd::SGAlignPath
end

struct SGAlignPaired <: SGAlignContainer
   fwd::SGAlignPath
   rev::SGAlignPath
end

# Store unique path's for gene/node but ignore
# SGAlignScore which is irrelevant at this stage
Base.isequal( a::SGAlignSingle, b::SGAlignSingle ) = a.fwd == b.fwd
Base.isequal( a::SGAlignPaired, b::SGAlignPaired ) = a.fwd == b.fwd && a.rev == b.rev

Base.:(==)( a::SGAlignSingle, b::SGAlignSingle ) = a.fwd == b.fwd
Base.:(==)( a::SGAlignPaired, b::SGAlignPaired ) = a.fwd == b.fwd && a.rev == b.rev

Base.hash( v::SGAlignSingle ) = hash( v.fwd )
Base.hash( v::SGAlignSingle, h::UInt64 ) = hash( v.fwd, h )
Base.hash( v::SGAlignPaired ) = hash( v.fwd, hash( v.rev ) )
Base.hash( v::SGAlignPaired, h::UInt64 ) = hash( v.fwd, hash( v.rev, h ) )

ispaired( cont::SGAlignSingle ) = false
ispaired( cont::SGAlignPaired ) = true

# check if one alignment path is a subset of another
function Base.in( first::SGAlignPath, last::SGAlignPath )
   for i in 1:(length(last)-length(first))+1
      offset = i - 1
      match = true
      for j in 1:length(first)
         if !(first[j] == last[offset+j])
            match = false
            break
         end
      end
      if match
         return true
      end
   end
   return false
end


# reimplementation of to edge code in events.jl
# we require not just the nodes to exist in the intset
# but also that they are adjacent in the intset
function Base.in( first, last, is::BitSet )
   if first in is &&
      last in is
      for i in first+1:last-1
         if i in is
            return false
         end
      end
      return true
   end
   return false
end

Base.in( edge::I, is::BitSet ) where I <: IntervalValue = in( edge.first, edge.last, is )

# This looks for an edge in a Vector of BitSets using ^
function Base.in( edge::I, viset::Vector{BitSet} ) where I <: IntervalValue
   for iset in viset
      if edge in iset
         return true
      end
   end
   false
end

# For events.jl:
function inall( edge::I, viset::Vector{BitSet} ) where I <: IntervalValue
   inone = false
   for iset in viset
      if edge.first >= first(iset) && edge.last <= last(iset)
         if !(edge in iset)
            return false
         end
         inone = true
      end
   end
   inone ? true : false
end

# each edge in the path must exist in the int set to be true
function Base.in( path::SGAlignPath, is::BitSet )
   if length(path) == 1
      if !(path[1].node in is)
         return false
      end
   else
      for i in 1:length(path)-1
         if !in( path[i].node, path[i+1].node, is )
            return false
         end
      end
   end
   return true
end

Base.in( aln::SGAlignSingle, is::BitSet ) = in( aln.fwd, is )
Base.in( aln::SGAlignPaired, is::BitSet ) = in( aln.fwd, is ) && in( aln.rev, is )

# This is where we count reads for nodes/edges/circular-edges/effective_lengths
# bias is an adjusting multiplier for nodecounts to bring them to the same level
# as junction counts, which should always be at a lower level
mutable struct SpliceGraphQuant{C <: SGAlignContainer, R <: ReadCounter}
   node::Vector{R}
   edge::IntervalMap{ExonInt,R}
   long::Dict{C,R}
   circ::Dict{Tuple{ExonInt,ExonInt},R}
   leng::Vector{Float64}
   bias::Float64

   function SpliceGraphQuant{C,R}( sg::SpliceGraph ) where {C <: SGAlignContainer, R <: ReadCounter}
      return new( R[zero(R) for x in 1:length(sg.nodelen)],
             IntervalMap{ExonInt,R}(),
             Dict{C,R}(),
             Dict{Tuple{ExonInt,ExonInt},R}(),
             zeros( length(sg.nodelen) ), 1.0 )
   end
end

abstract type EquivalenceClass end

const ClassInt  = Int32
const ClassType = Vector{ClassInt}

# type-stable isoform compatibility class
mutable struct IsoCompat <: EquivalenceClass
   class::ClassType
   prop::Vector{Float64}
   prop_sum::Float64
   count::Float64
   isdone::Bool

   function IsoCompat( arr::Vector{Int}, count::Float64 )
      return new( arr, ones( length(arr) ) / length(arr), 1.0, count, false )
   end
end

# Here we store whole graphome quantification
struct GraphLibQuant{C <: SGAlignContainer, R <: ReadCounter}
   tpm::Vector{Float64}                 # isoform tpm
   count::Vector{Float64}               # isoform counts
   length::Vector{Int64}                # isoform lengths
   geneidx::Vector{Int64}               # 0-based offset of isoforms for gene
   genetpm::Vector{Float64}             # sum of isoform tpms for each gene
   quant::Vector{SpliceGraphQuant{C,R}} # splice graph quant structs for gene
   classes::Vector{IsoCompat}           # isoform compatibility classes

   function GraphLibQuant{C,R}( lib::GraphLib ) where {C <: SGAlignContainer, R <: ReadCounter}
      isonum = zeros(Int, length(lib.graphs))
      for i in 1:length(lib.graphs)
         isonum[i] = length(lib.graphs[i].annoname)
      end
      isolen  = sum( isonum )
      tpm     = zeros( isolen )
      count   = zeros( isolen )
      leng    = zeros( isolen )
      geneoff = zeros(Int, length(lib.graphs))
      genetpm = zeros( length(lib.graphs) )
      quant   = Vector{SpliceGraphQuant}(undef, length(lib.graphs))
      classes = Vector{IsoCompat}()
      cumul_i = 0
      for i in 1:length(lib.graphs)
         geneoff[i] = cumul_i
         for j in 1:length(lib.graphs[i].annopath)
            leng[cumul_i+j] = sum(map( y->lib.graphs[i].nodelen[y], collect(lib.graphs[i].annopath[j]) ))
         end
         cumul_i += isonum[i]
         quant[i] = SpliceGraphQuant{C,R}( lib.graphs[i] )
      end
      new( tpm, count, leng, geneoff, genetpm, quant, classes )
   end
end

# Adjust each MultiCompat raw value and assign its adjusted to .count
function adjust!( sgquant::SpliceGraphQuant{C,R}, mod::B, func=adjust! ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel}
   for node in sgquant.node
      func( node, mod )
   end
   for edg in sgquant.edge
      func( edg.value, mod )
   end
   for lng in values(sgquant.long)
      func( lng, mod )
   end
   for cir in values(sgquant.circ)
      func( cir, mod )
   end
   sgquant
end

primer_adjust!( sgquant::SpliceGraphQuant{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( sgquant, mod, primer_adjust! )
gc_adjust!( sgquant::SpliceGraphQuant{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( sgquant, mod, gc_adjust! )

# Adjust each sgquant storage value
function adjust!( gquant::GraphLibQuant{C,R}, mod::B, func=adjust! ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel}
   for sgq in gquant.quant
      func( sgq, mod )
   end
   gquant
end

primer_adjust!( gquant::GraphLibQuant{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( gquant, mod, primer_adjust! )
gc_adjust!( gquant::GraphLibQuant{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( gquant, mod, gc_adjust! )

function gc_normalize!( mod::JointBiasMod, lib::GraphLib, quant::GraphLibQuant )
   for i in 1:length(lib.graphs)
      sg  = lib.graphs[i]
      idx = quant.geneidx[i]
      for j in 1:length(sg.annobias)
         add_to_back!( mod.gcback, sg.annobias[j], quant.tpm[j+idx] )
      end
   end
   foresum = sum(mod.gcfore)
   backsum = sum(mod.gcback)
   @inbounds for i in 1:length(mod.gcfore)
      @fastmath mod.gcfore[i] /= foresum
      @fastmath mod.gcback[i] /= backsum
   end
   mod
end

@inbounds function set_gene_tpm!( quant::GraphLibQuant, lib::GraphLib )
   for i in 1:length(lib.graphs)
      tpm = 0.0
      curidx = quant.geneidx[i]
      for j in 1:length( lib.graphs[i].annoname )
         tpm += quant.tpm[ curidx + j ]
      end
      quant.genetpm[i] = tpm
   end
end

assign_path!( graphq::GraphLibQuant{SGAlignSingle}, path::SGAlignSingle,
              value, temp_iset::BitSet ) = assign_path!( graphq, path, value )

# This function re-count!'s SGAlignPaths which were not count!'ed originally
@inbounds function assign_path!( graphq::GraphLibQuant{SGAlignSingle}, path::SGAlignSingle, value )
   init_gene = path.fwd[1].gene
   sgquant = graphq.quant[ init_gene ]

   if length(path.fwd) == 1
      push!( sgquant.node[ path.fwd[1].node ], value )
   else
      for i in 1:length(path.fwd)-1
         lnode = path.fwd[i].node
         rnode = path.fwd[i+1].node
         if lnode < rnode
            interv = IntervalTrees.Interval{ExonInt}( lnode, rnode )
            pushzero!( sgquant.edge, interv, value )
         elseif lnode >= rnode
            pushzero!( sgquant.circ, (lnode,rnode), value )
         end
      end
   end
end

# This function re-counts paired SGAlignPaths that weren't counted because they were too long, disjointed, or both
@inbounds function assign_path!( graphq::GraphLibQuant{SGAlignPaired}, path::SGAlignPaired, value, temp_iset::BitSet=BitSet() )
   init_gene = path.fwd[1].gene
   sgquant = graphq.quant[ init_gene ]

   if length(path.fwd) == 1
      push!( sgquant.node[ path.fwd[1].node ], value )
      push!( temp_iset, path.fwd[1].node )
   else
      for i in 1:length(path.fwd)-1
         lnode = path.fwd[i].node
         rnode = path.fwd[i+1].node
         push!(temp_iset, lnode)
         push!(temp_iset, rnode)
         if lnode < rnode
            interv = IntervalTrees.Interval{ExonInt}( lnode, rnode )
            pushzero!( sgquant.edge, interv, value )
         elseif lnode >= rnode
            pushzero!( sgquant.circ, (lnode,rnode), value )
         end
      end
   end

   if length(path.rev) == 1
      if !(path.rev[1].node in temp_iset)
         push!( sgquant.node[ path.rev[1].node ], value )
      end
   else
      for i in 1:length(path.rev)-1
         lnode = path.rev[i].node
         rnode = path.rev[i+1].node
         (lnode in temp_iset && rnode in temp_iset) && continue
         if lnode < rnode
            interv = IntervalTrees.Interval{ExonInt}( lnode, rnode )
            pushzero!( sgquant.edge, interv, value )
         elseif lnode >= rnode
            pushzero!( sgquant.circ, (lnode,rnode), value )
         end
      end
   end
   empty!(temp_iset)
end

# build compatibility classes from SpliceGraphQuant node and edge "counts"
@inbounds function build_equivalence_classes!( graphq::GraphLibQuant, lib::GraphLib; assign_long=true )
   # initialize re-used data structures
   iset = BitSet()
   # resize!(iset.bits, 64)
   temp = Dict{Vector{Int},Float64}()

   @inline function handle_iset_value!( idx, value )
      if length(iset) == 1
         graphq.count[first(iset) + idx] += value
      elseif length(iset) > 1
         arr = collect(iset)
         if haskey(temp, arr)
            temp[arr] += value
         else
            temp[arr] = value
         end
      end
   end

   # go through splice graphs (to re-use temp data structures)
   for i in 1:length(lib.graphs)
      # go through single node mapping reads
      for n in 1:length(graphq.quant[i].node)
         for p in 1:length(lib.graphs[i].annopath)
            if n in lib.graphs[i].annopath[p]
               push!(iset, p)
            end
         end
         handle_iset_value!( graphq.geneidx[i], get(graphq.quant[i].node[n]) )
         empty!(iset)
      end
      # go through edges
      for edg in graphq.quant[i].edge
         for p in 1:length(lib.graphs[i].annopath)
            if edg in lib.graphs[i].annopath[p]
               push!(iset, p)
            end
         end
         handle_iset_value!( graphq.geneidx[i], get(edg.value) )
         empty!(iset)
      end
      # go through multi-edge crossing paths
      for align in keys(graphq.quant[i].long)
         for p in 1:length(lib.graphs[i].annopath)
            if align in lib.graphs[i].annopath[p]
               push!(iset, p)
            end
         end
         handle_iset_value!( graphq.geneidx[i], get(graphq.quant[i].long[align]) )
         empty!(iset)
         if assign_long
            assign_path!( graphq, align, get(graphq.quant[i].long[align]), iset )
         end
      end
      # create compatibility classes for isoforms
      for arr in keys(temp)
         if temp[arr] > 0.0
            push!( graphq.classes, IsoCompat(arr .+ graphq.geneidx[i], temp[arr]) )
         end
      end
      empty!(temp)
   end
end


# This function takes a vector of alignments and returns the compatibility class
# and a boolean describing whether or not we can assign the read through EM
@inbounds function equivalence_class( graphq::GraphLibQuant,
                                      lib::GraphLib,
                                      cont::Vector{C},
                                      temp_iset::BitSet=BitSet() ) where C <: SGAlignContainer
   matches_all = true
   for path in cont
      g = path.fwd[1].gene
      has_match = false
      for p in 1:length(lib.graphs[g].annopath)
         if path in lib.graphs[g].annopath[p]
            push!(temp_iset, p + graphq.geneidx[g])
            has_match = true
         end
      end
      if !has_match
         matches_all = false
         empty!(temp_iset)
         return ClassType(),false
      end
   end
   retval = collect(ClassInt, temp_iset)
   empty!(temp_iset)
   retval,matches_all
end


## MultiAlignment containers & Equivalence classes

mutable struct MultiCompat{R} <: EquivalenceClass
   class::ClassType
   prop::Vector{Float64}
   prop_sum::Float64
   count::Float64
   isdone::Bool
   postassign::Bool
   raw::R

   function MultiCompat{R}( graphq::GraphLibQuant{C,R}, lib::GraphLib, cont::Vector{C}, temp_iset::BitSet=BitSet() ) where {C <: SGAlignContainer, R <: ReadCounter}
      class,matches = equivalence_class( graphq, lib, cont, temp_iset )
      new( class, ones(length(class)) / length(class), 1.0, 1.0, !matches, !matches, zero(R) )
   end
end

# This hash structure stores multi-mapping
# equivalence classes
struct MultiMapping{C <: SGAlignContainer, R <: ReadCounter}
   map::Dict{Vector{C},MultiCompat{R}}
   iset::BitSet

   function MultiMapping{C,R}() where {C <: SGAlignContainer, R <: ReadCounter}
      new( Dict{Vector{C},MultiCompat{R}}(), BitSet() )
   end
end

function total_multi( mm::MultiMapping )
   total = 0.0
   for v in values(mm.map)
      total += v.count
   end
   total
end

function Base.push!( ambig::MultiMapping{SGAlignSingle,R}, alns::Vector{SGAlignment}, value,
                     graphq::GraphLibQuant, lib::GraphLib ) where R <: ReadCounter
   cont = Vector{SGAlignSingle}(undef, length(alns))
   for i in 1:length(alns)
      cont[i] = SGAlignSingle( alns[i].path )
   end
   if haskey( ambig.map, cont )
      push!( ambig.map[cont].raw, value )
   else
      mc = MultiCompat{R}( graphq, lib, cont, ambig.iset )
      push!( mc.raw, value )
      ambig.map[cont] = mc
   end
end

function Base.push!( ambig::MultiMapping{SGAlignPaired,R}, fwd::Vector{SGAlignment}, rev::Vector{SGAlignment}, value,
                     graphq::GraphLibQuant, lib::GraphLib ) where R <: ReadCounter
   cont = Vector{SGAlignPaired}( undef, length(fwd) )
   for i in 1:length(fwd)
      cont[i] = SGAlignPaired( fwd[i].path, rev[i].path )
   end
   if haskey( ambig.map, cont )
      push!( ambig.map[cont].raw, value )
   else
      mc = MultiCompat{R}( graphq, lib, cont, ambig.iset )
      push!( mc.raw, value )
      ambig.map[cont] = mc
   end
end

# Adjust each MultiCompat raw value and assign its adjusted to .count
function adjust!( multi::MultiMapping{C,R}, mod::B, func=adjust! ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel}
   for mc in values(multi.map)
      func( mc.raw, mod )
      mc.count = get(mc.raw)
   end
   multi
end

primer_adjust!( multi::MultiMapping{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( multi, mod, primer_adjust! )
gc_adjust!( multi::MultiMapping{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( multi, mod, gc_adjust! )

# for a given alignment container, set the temp_iset in multi to the compatibility class
@inbounds function set_equivalence_class!( multi::MultiMapping{C,R}, graphq::GraphLibQuant,
                                           lib::GraphLib, aln::C ) where {C <: SGAlignContainer, R <: ReadCounter}
   empty!(multi.iset)
   g = aln.fwd[1].gene
   for p in 1:length(lib.graphs[g].annopath)
      if aln in lib.graphs[g].annopath[p]
         push!(multi.iset, p + graphq.geneidx[g])
      end
   end
end

function get_class_value( iset::BitSet, compat::MultiCompat )
   value = 0.0
   @inbounds for i in 1:length(compat.class)
      if compat.class[i] in iset
         value += compat.prop[i]
      end
   end
   value
end

function assign_ambig!( graphq::GraphLibQuant{C,R}, lib::GraphLib, ambig::MultiMapping{C,R} ) where {C <: SGAlignContainer, R <: ReadCounter}
   @inbounds for (vp,mc) in ambig.map
      # store relative weight of each alignment path
      vals = zeros(length(vp))

      if mc.postassign  # not all paths match annotated transcripts, use gene level values to assign
         total_tpm = 0.0
         for i in 1:length(vp)
            vals[i] = graphq.genetpm[ vp[i].fwd[1].gene ]
            total_tpm += vals[i]
         end
         @fastmath vals /= sum(vals)
      else # assign based on transcripts
         # since two alignments for a read can match the same transcript
         # lets sum the prop values of the MultiCompat for each alignment
         # then we normalize by the total sum of those, and assign reads accordingly.
         # Example: imagine a single read matches the same transcript in two positions
         # that creates a single MultiCompat object that has two alignments hitting the same
         # transcript class.  This method would assign that prop to both, and then re-normalize
         # so that the total count of the MultiCompat is assigned for the sum of the alignments
         vals = zeros(length(vp))
         for i in 1:length(vp)
            set_equivalence_class!( ambig, graphq, lib, vp[i] )
            vals[i] = get_class_value( ambig.iset, mc )
         end
         @fastmath vals /= sum(vals)
      end

      # now assign alignment paths proportional to their weight
      for i in 1:length(vp)
         assign_path!( graphq, vp[i], vals[i] * mc.count, ambig.iset )
      end
   end
end

# Count a single alignment in quant
function count!( graphq::GraphLibQuant{C,R}, align::SGAlignment, value ) where {C <: SGAlignContainer, R <: ReadCounter}
   align.isvalid == true || return
   init_gene = align.path[1].gene
   sgquant   = graphq.quant[ init_gene ]

   # single node support
   if length(align.path) == 1
      # access node -> [ SGNode( gene, *node* ) ]
      push!( sgquant.node[ align.path[1].node ], value )
   else # multi-node mapping read
      # if any node in path maps to other gene, ignore read
      for n in 1:length(align.path)
         if align.path[n].gene != init_gene
            return
         end
      end

      # add a single edge support
      if length(align.path) == 2
         lnode = align.path[1].node
         rnode = align.path[2].node
         if lnode < rnode
            interv = IntervalTrees.Interval{ExonInt}( lnode, rnode )
            pushzero!( sgquant.edge, interv, value )
         elseif lnode >= rnode
            pushzero!( sgquant.circ, (lnode,rnode), value )
         end
      else # or we have long path, add long count
         curkey = SGAlignSingle(align.path)
         if haskey(sgquant.long, curkey)
            push!( sgquant.long[curkey], value )
         else
            entry = zero(R)
            push!( entry, value )
            sgquant.long[curkey] = entry
         end
      end
   end
end

## Extension of count! for paired end counting
function count!( graphq::GraphLibQuant{C,R}, fwd::SGAlignment, rev::SGAlignment, value ) where {C <: SGAlignContainer, R <: ReadCounter}
   (fwd.isvalid == true && rev.isvalid == true) || return
   init_gene = fwd.path[1].gene
   rev_gene  = rev.path[1].gene
   (init_gene == rev_gene) || return
   sgquant = graphq.quant[ init_gene ]
   singlebool = false

   if length(fwd.path) <= length(rev.path)
      if fwd.path in rev.path
         singlebool = true
         single = rev.path
      end
   else
      if rev.path in fwd.path
         singlebool = true
         single = fwd.path
      end
   end

   # single node support
   if singlebool
      if length(single) == 1
         # access node -> [ SGNode( gene, *node* ) ]
         push!( sgquant.node[ single[1].node ], value )
      else # multi-node mapping read
         # if any node in path maps to other gene, ignore read
         for n in 1:length(single)
            if single[n].gene != init_gene
               return
            end
         end

         # add a single edge support
         if length(single) == 2
            lnode = single[1].node
            rnode = single[2].node
            if lnode < rnode
               interv = IntervalTrees.Interval{ExonInt}( lnode, rnode )
               pushzero!( sgquant.edge, interv, value )
            elseif lnode >= rnode
               pushzero!( sgquant.circ, (lnode,rnode), value )
            end
         end
      end
   else
       # or we have two paths, add long count
       curkey = SGAlignPaired(fwd.path, rev.path)
       if haskey(sgquant.long, curkey)
          push!( sgquant.long[curkey], value )
       else
          entry = zero(R)
          push!( entry, value )
          sgquant.long[curkey] = entry
       end
   end
end

@inline function calculate_tpm!( quant::GraphLibQuant, counts::Vector{Float64}=quant.count; readlen::Int64=50, sig::Int64=1 )
   for i in 1:length(counts)
      @fastmath quant.tpm[ i ] = counts[i] / max( (quant.length[i] - readlen), 1.0 )
   end
   rpk_sum = max( sum( quant.tpm ), 1.0 )
   for i in 1:length(quant.tpm)
      if sig > 0
         @fastmath quant.tpm[i] = round( quant.tpm[i] * SCALING_FACTOR / rpk_sum, digits=sig )
      else
         @fastmath quant.tpm[i] = ( quant.tpm[i] * SCALING_FACTOR / rpk_sum )
      end
   end
end

calculate_bias!( sgquant::SpliceGraphQuant ) = 1.0
global_bias( graphq::GraphLibQuant ) = 1.0,0.0

# Every exon-exon junction has an effective mappable space of readlength - K*2.
# Since we only count nodes that don't map over edges, we can give
# give each junction an effective_length of one, and now the space
# inside of a node for which a read could map without overlapping a junction
# is put into units of junction derived effective_length
# kadj is minimum of (readlength - minalignlen) and k-1
@inline function eff_length( node, sg::SpliceGraph, eff_len::Int, kadj::Int )
   len = sg.nodelen[node] + (istxstart( sg.edgetype[node] ) ? 0 : kadj) +
                            (istxstop( sg.edgetype[node+1] ) ? 0 : kadj)
   @fastmath len / eff_len
end

function eff_lengths!( sg::SpliceGraph, sgquant::SpliceGraphQuant, eff_len::Int, kadj::Int )
   for i in 1:length( sg.nodelen )
      sgquant.leng[i] = sg.nodelen[i] #eff_length( i, sg, eff_len, kadj )
   end
end

function effective_lengths!( lib::GraphLib, graphq::GraphLibQuant, eff_len::Int, kadj::Int )
   for i in 1:length( lib.graphs )
      eff_lengths!( lib.graphs[i], graphq.quant[i], eff_len, kadj )
   end
end

function unsafe_copy!( dest::Vector{T}, src::Vector{T}; indx_shift=0 ) where T <: Number
   @inbounds for i in 1:length(src)
      dest[i+indx_shift] = src[i]
   end
end

function copy_isdone!( dest::Vector{T}, src::Vector{T}, denominator::Float64=1.0, sig::Int=4 ) where T <: AbstractFloat
   isdifferent = false
   for i in 1:length(dest)
      val = round( src[i] / denominator, digits=sig )
      if val != dest[i]
         isdifferent = true
      end
      dest[i] = isnan(val) ? 0.0 : val
   end
   !isdifferent || sum(denominator) == 0
end

# This function performs expectation maximization
# and sets the quant.tpm array as the proposed expression set
# at each iteration.   It also requires that calculate_tpm! be
# run initially once prior to gene_em!() call.

function gene_em!( quant::GraphLibQuant, ambig::MultiMapping{C,R};
                   it::Int64=1, maxit::Int64=1000,
                   sig::Int64=4, readlen::Int64=50 ) where {C <: SGAlignContainer, R <: ReadCounter}

   count_temp = ones(length(quant.count))
   tpm_temp   = ones(length(quant.count))
   prop_temp  = zeros( 500 )
   uniqsum    = sum(quant.count)
   ambigsum   = length(ambig.map)

   @inline function maximize_and_assign!( eq::E, maximize::Bool=true ) where E <: EquivalenceClass
      eq.isdone && return
      if maximize
         eq.prop_sum = 0.0
         c = 1
         for i in eq.class
            cur_tpm = quant.tpm[i]
            prop_temp[c] = cur_tpm
            eq.prop_sum += cur_tpm
            c += 1
         end
         eq.isdone = copy_isdone!( eq.prop, prop_temp, eq.prop_sum )
      end
      c = 1
      for i in eq.class
         if eq.isdone
            quant.count[i] += eq.prop[c] * eq.count
         end
         count_temp[i] += eq.prop[c] * eq.count
         c += 1
      end
   end

   while tpm_temp != quant.tpm && it < maxit

      unsafe_copy!( count_temp, quant.count )
      unsafe_copy!( tpm_temp, quant.tpm )

      for eq in values(ambig.map)
         maximize_and_assign!( eq, (it > 1) )
      end
      for eq in quant.classes
         maximize_and_assign!( eq, (it > 1) )
      end

      calculate_tpm!( quant, count_temp, sig=sig, readlen=readlen )

      it += 1
   end
   it
end


function output_tpm( file, lib::GraphLib, gquant::GraphLibQuant )
   geneio = open( file * ".gene.tpm.gz", "w" )
   genestream = ZlibDeflateOutputStream( geneio )
   isoio = open( file * ".isoform.tpm.gz", "w" )
   isostream = ZlibDeflateOutputStream( isoio )
   output_tpm_header( genestream, "Gene" )
   output_tpm_header( isostream, "Isoform" )
   for i in 1:length(lib.names)
      idx = gquant.geneidx[i]
      genetpm = 0.0
      genecnt = 0.0
      for j in 1:length(lib.graphs[i].annoname)
         genetpm += gquant.tpm[j+idx]
         genecnt += gquant.count[j+idx]
         tab_write( isostream, lib.graphs[i].annoname[j] )
         tab_write( isostream, string(gquant.tpm[j+idx]) )
         tab_write( isostream, string(gquant.count[j+idx]) )
         write( isostream, '\n' )
      end
      tab_write( genestream, lib.names[i] )
      tab_write( genestream, string(genetpm) )
      write( genestream, string(genecnt) )
      write( genestream, '\n' )
   end
   close( isostream )
   close( isoio )
   close( genestream )
   close( geneio )
end
