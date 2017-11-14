
const SCALING_FACTOR = 1_000_000

# use specialty type for "read counts"
# with basic function: 'sum' that can be
# overridden for more complex count bias models
const ReadCount = Float64

# This is where we count reads for nodes/edges/circular-edges/effective_lengths
# bias is an adjusting multiplier for nodecounts to bring them to the same level
# as junction counts, which should always be at a lower level, e.g. bias < 1
mutable struct SpliceGraphQuant
   node::Vector{ReadCount}
   edge::IntervalMap{ExonInt,ReadCount}
   long::Dict{SGAlignPath,ReadCount}
   circ::Dict{Tuple{ExonInt,ExonInt},ReadCount}
   leng::Vector{Float64}
   bias::Float64
end

# Default constructer
SpliceGraphQuant() = SpliceGraphQuant( Vector{ReadCount}(),
                                       IntervalMap{ExonInt,ReadCount}(),
                                       Dict{SGAlignPath,ReadCount}(),
                                       Dict{Tuple{ExonInt,ExonInt},ReadCount}(),
                                       Vector{Float64}(), 1.0 )

SpliceGraphQuant( sg::SpliceGraph ) = SpliceGraphQuant( zeros(ReadCount, length(sg.nodelen) ),
                                                        IntervalMap{ExonInt,ReadCount}(),
                                                        Dict{SGAlignPath,ReadCount}(),
                                                        Dict{Tuple{ExonInt,ExonInt},ReadCount}(),
                                                        zeros( length(sg.nodelen) ), 1.0 )

abstract type EquivalenceClass end

# type-stable isoform compatibility class
mutable struct IsoCompat <: EquivalenceClass
   class::Vector{Int16}
   prop::Vector{Float64}
   prop_sum::Float64
   count::Float64
   isdone::Bool

   function IsoCompat( arr::Vector{Int}, count::Float64 )
      return IsoCompat( arr, ones( length(arr) ) / length(arr), 1.0, count, false )
   end
end


IsoQuant( sg::SpliceGraph ) = IsoQuant( zeros( length(sg.annoname) ),
                                        zeros( length(sg.annoname) ),
                                        map( x->Int32(sum(map( y->sg.nodelen[y], collect(x) ))), sg.annopath ),
                                        Vector{IsoCompat}(), true )

function Base.in( path::SGAlignPath, is::IntSet )
   for i in 1:length(path)
      if !(path[i].node in is)
         return false
      end
   end
   return true
end

# build isoquant compatibility classes from SpliceGraphQuant node and edge "counts"
@inbounds function fill_isoquant!( graphq::GraphLibQuant, lib::GraphLib; assign_edge=true )
   # initialize re-used data structures
   iset = IntSet()
   resize!(iset.bits, 64)
   temp = Dict{Tuple,Float64}()

   @inline function handle_iset_value!( value )
      if length(iset) == 1
         graphq.compat[i].count[first(iset)] += sum(value)
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
         handle_iset_value!( sum(graphq.quant[i].node[n]) )
         empty!(iset)
      end
      # go through edges
      for edg in graphq.quant[i].edge
         for p in 1:length(lib.graphs[i].annopath)
            if edg.first in lib.graphs[i].annopath[p] &&
               edg.last  in lib.graphs[i].annopath[p]
               push!(iset, p)
            end
         end
         handle_iset_value!( sum(edg.value) )
         empty!(iset)
      end
      # go through multi-edge crossing paths
      for align in keys(graphq.quant[i].long)
         for p in 1:length(lib.graphs[i].annopath)
            if align in lib.graphs[i].annopath[p]
               push!(iset, p)
            end
         end
         handle_iset_value!( sum(graphq.quant[i].long[align]) )
         empty!(iset)
         if assign_edge
            #TODO!!!
         end
      end
      # create compatibility classes for isoforms
      for arr in keys(temp) #TODO
         push!( graphq.classes, IsoCompat(arr, temp[arr]) )
      end
      empty!(temp)
   end
end

# Here we store whole graphome quantification
struct GraphLibQuant
   tpm::Vector{Float64}
   count::Vector{Float64}
   length::Vector{Int64}
   geneidx::Vector{Int64}
   quant::Vector{SpliceGraphQuant}
   classes::Vector{IsoCompat}
end

function GraphLibQuant( lib::GraphLib )
   isonum  = map( x->length(x.annopath), lib.graphs )
   tpm     = zeros( isolen )
   count   = zeros( isolen )
   length  = zeros( isolen )
   geneoff = zeros(Int, length(lib.graphs))
   quant   = Vector{SpliceGraphQuant}( length(lib.graphs) )
   cumul_i = 0
   for i in 1:length(lib.graphs)
      geneoff[i] = cumul_i
      cumul_i += isonum[i]
      # START HERE!
      quant[i] = SpliceGraphQuant( lib.graphs[i] )
   end
   GraphLibQuant( tpm, isoq, quant )
end

#=
@inline function calculate_tpm!( quant::GraphLibQuant, counts::Vector{Float64}=quant.count; readlen::Int64=50, sig::Int64=1 )
   for i in 1:length(counts)
      @fastmath quant.tpm[ i ] = counts[i] / max( (quant.length[i] - readlen), 1.0 )
   end
   const rpk_sum = max( sum( quant.tpm ), 1.0 )
   for i in 1:length(quant.tpm)
      if sig > 0
         @fastmath quant.tpm[i] = round( quant.tpm[i] * SCALING_FACTOR / rpk_sum, sig )
      else
         @fastmath quant.tpm[i] = ( quant.tpm[i] * SCALING_FACTOR / rpk_sum )
      end
   end
end
=#

@inbounds function calculate_tpm!( quant::GraphLibQuant; readlen::Int64=50, sig::Int64=1 )
   for i in 1:length(quant.compat)
      quant.tpm[i] = 0.0
      for j in 1:length(quant.compat[i].rpk)
         @fastmath quant.tpm[i] += quant.compat[i].rpk[j] / max( quant.compat[i].length[j] - readlen, 1.0 )
      end
   end
   const rpk_sum = max( sum( quant.tpm ), 1.0 )
   for i in 1:length(quant.tpm)
      if sig > 0
         @fastmath quant.tpm[i] = round( quant.tpm[i] * SCALING_FACTOR / rpk_sum, sig )
      else
         @fastmath quant.tpm[i] = ( quant.tpm[i] * SCALING_FACTOR / rpk_sum )
      end
   end
end

const MultiAln = Vector{SGAlignPath}

mutable struct MultiCompat <: EquivalenceClass
   class::Vector{Int16}
   prop::Vector{Float64}
   prop_sum::Float64
   count::Float64
   raw::ReadCount
   isdone::Bool
end

function MultiAssign( 

# This hash structure stores multi-mapping
# equivalence classes
const MultiCompat = Dict{MultiAln,MultiAssign}

function Base.push!( ambig::MultiCompat, aligns::Vector{SGAlignment} )
   
end

#=
mutable struct Multimap
   align::Vector{SGAlignment}
   prop::Vector{Float64}
   prop_sum::Float64
end


Multimap( aligns::Vector{SGAlignment} ) = length(aligns) >= 1 ? 
                                    Multimap( aligns, ones(length(aligns)) / length(aligns), 1.0 ) :
                                    Multimap( aligns, Float64[], 0.0 )
=#

function assign_ambig!( graphq::GraphLibQuant, ambig::Vector{Multimap}; ispaired::Bool=false )
   for mm in ambig
      i = 1
      while i <= length(mm.prop)
         if ispaired
            (i == length(mm.prop)) && break
            count!( graphq, mm.align[i], mm.align[i+1], val=mm.prop[i] )
            i += 1
         else
            count!( graphq, mm.align[i], val=mm.prop[i] )
         end
         i += 1
      end
   end
end

function count!( graphq::GraphLibQuant, align::SGAlignment; val::Float64=1.0 )
   align.isvalid == true || return
   init_gene = align.path[1].gene
   sgquant   = graphq.quant[ init_gene ]

   # single node support
   if length(align.path) == 1
      # access node -> [ SGNode( gene, *node* ) ]
      sgquant.node[ align.path[1].node ] += val
   else # multi-node mapping read
      # if any node in path maps to other gene, ignore read
      for n in 1:length(align.path)
         if align.path[n].gene != init_gene
            return
         end
      end

      # add a single edge support
      if length(align.path) == 2
         const lnode = align.path[1].node
         const rnode = align.path[2].node
         if lnode < rnode
            interv = Interval{ExonInt}( lnode, rnode )
            sgquant.edge[ interv ] = get( sgquant.edge, interv, IntervalValue(0,0,0.0) ).value + val
         elseif lnode >= rnode
            sgquant.circ[ (lnode, rnode) ] = get( sgquant.circ, (lnode,rnode), 0.0) + val
         end         
      else # or we have long path, add long count
         if haskey(sgquant.long, align.path)
            sgquant.long[align.path] += val
         else
            sgquant.long[align.path] = val
         end
      end
   end
end

function calculate_bias!( sgquant::SpliceGraphQuant )
   edgecnt = 0.0
   for edgev in sgquant.edge
      edgecnt += edgev.value
   end
   @fastmath edgelevel = edgecnt / length(sgquant.edge)
   nodecnt = 0.0
   for nodev in sgquant.node
      nodecnt += nodev
   end
   @fastmath nodelevel = nodecnt / sum( sgquant.leng )
   # never down-weight junction reads, only exon-body reads
   bias = @fastmath min( edgelevel / nodelevel, 1.0 )
   sgquant.bias = bias
   bias
end

function global_bias( graphq::GraphLibQuant )
   bias_ave = 0.0
   bias_var = 0.0
   n = 1
   for sgq in graphq.quant
      curbias = calculate_bias!( sgq )
      if curbias > 10
         println( sgq )
      end
      old = bias_ave
      @fastmath bias_ave += (curbias - bias_ave) / n
      if n > 1
        @fastmath bias_var += (curbias - old)*(curbias - bias_ave)
      end
      n += 1
   end
   bias_ave,bias_var
end

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
      sgquant.leng[i] = eff_length( i, sg, eff_len, kadj )
   end
end

function effective_lengths!( lib::GraphLib, graphq::GraphLibQuant, eff_len::Int, kadj::Int )
   for i in 1:length( lib.graphs )
      eff_lengths!( lib.graphs[i], graphq.quant[i], eff_len, kadj )
   end
end

function Base.unsafe_copy!{T <: Number}( dest::Vector{T}, src::Vector{T}; indx_shift=0 )
   for i in 1:length(src)
      dest[i+indx_shift] = src[i]
   end
end

# This function performs expectation maximization
# and sets the quant.tpm array as the proposed expression set
# at each iteration.   It also requires that calculate_tpm! be 
# run initially once prior to gene_em!() call.

function gene_em!( quant::GraphLibQuant, ambig::Vector{Multimap};
                   it::Int64=1, maxit::Int64=1000, sig::Int64=1, readlen::Int64=50 )

   const count_temp = ones(length(quant.count))
   const tpm_temp   = ones(length(quant.count))
   const uniqsum    = sum(quant.count)
   const ambigsum   = length(ambig)

   while tpm_temp != quant.tpm && it < maxit

      unsafe_copy!( count_temp, quant.count )
      unsafe_copy!( tpm_temp, quant.tpm )

      for mm in ambig
         if it > 1 # Maximization
            mm.prop_sum = 0.0
            for ai in 1:length(mm.prop)
               const init_gene = mm.align[ai].path[1].gene
               const init_tpm  = quant.tpm[ init_gene ] * max( quant.length[ init_gene ] - readlen, 1.0 )
               mm.prop[ai] = init_tpm
               @fastmath mm.prop_sum += init_tpm
            end
         end

         for ai in 1:length(mm.align)
            const init_gene = mm.align[ai].path[1].gene
            @fastmath const prop = mm.prop[ai] / mm.prop_sum
            mm.prop[ai] = isnan(prop) ? 0.0 : prop
            @fastmath count_temp[ init_gene ] += mm.prop[ai]
         end
      end

      calculate_tpm!( quant, count_temp, sig=sig, readlen=readlen ) # Expectation
 
      it += 1
   end
   it
end

function gene_em!( quant::GraphLibQuant, ambig::MultiCompat;
                   it::Int64=1, maxit::Int64=1000, sig::Int64=1, readlen::Int64=50 )

   const tpm_temp   = ones(length(quant.count))
   const uniqsum    = sum(quant.count)
   const ambigsum   = length(ambig)

   while tpm_temp != quant.tpm && it < maxit

      unsafe_copy!( tpm_temp, quant.tpm )

      for (k,v) in ambig
         if it > 1 # Maximization
            v.prop_sum = 0.0
            for ai in 1:length(k)
               const init_gene = mm.align[ai].path[1].gene
               const init_tpm  = quant.tpm[ init_gene ] * max( quant.length[ init_gene ] - readlen, 1.0 )
               mm.prop[ai] = init_tpm
               @fastmath mm.prop_sum += init_tpm
            end
         end

         for ai in 1:length(mm.align)
            const init_gene = mm.align[ai].path[1].gene
            @fastmath const prop = mm.prop[ai] / mm.prop_sum
            mm.prop[ai] = isnan(prop) ? 0.0 : prop
            @fastmath count_temp[ init_gene ] += mm.prop[ai]
         end
      end

      calculate_tpm!( quant, count_temp, sig=sig, readlen=readlen ) # Expectation
 
      it += 1
   end
   it
end


function output_tpm( file, lib::GraphLib, gquant::GraphLibQuant )
   io = open( file, "w" )
   stream = ZlibDeflateOutputStream( io )
   output_tpm_header( stream )
   for i in 1:length(lib.names)
      tab_write( stream, lib.names[i] )
      tab_write( stream, string(gquant.tpm[i]) )
      tab_write( stream, string(gquant.count[i]) )
      write( stream, '\n' )
   end
   close( stream )
   close( io )
end
