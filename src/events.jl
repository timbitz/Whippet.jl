
bitstype 8 EdgeMotif

Base.convert( ::Type{EdgeMotif}, motif::UInt8 ) = box(EdgeMotif, unbox(UInt8, motif ))
Base.convert( ::Type{UInt8}, motif::EdgeMotif ) = box(UInt8, unbox(EdgeMotif, motif ))

# OBLIGATE MOTIFS (provided the cognate nodes both have reads):

const TXST_MOTIF = convert(EdgeMotif, 0b000)
const TXEN_MOTIF = convert(EdgeMotif, 0b001)

const ALTF_MOTIF = convert(EdgeMotif, 0b010)
const ALTL_MOTIF = convert(EdgeMotif, 0b011)

# REQUIRES SPANNING EDGE:

const RETI_MOTIF = convert(EdgeMotif, 0b100) 
const SKIP_MOTIF = convert(EdgeMotif, 0b101)

const ALTD_MOTIF = convert(EdgeMotif, 0b110)
const ALTA_MOTIF = convert(EdgeMotif, 0b111)

# NO MOTIF
const NONE_MOTIF = convert(EdgeMotif, 0b1000)

const MOTIF_STRING = fill( "", 9 )
      MOTIF_STRING[ UInt8(TXST_MOTIF) + 1 ] = "TS"
      MOTIF_STRING[ UInt8(TXEN_MOTIF) + 1 ] = "TE"
      MOTIF_STRING[ UInt8(ALTF_MOTIF) + 1 ] = "AF"
      MOTIF_STRING[ UInt8(ALTL_MOTIF) + 1 ] = "AL"
      MOTIF_STRING[ UInt8(RETI_MOTIF) + 1 ] = "RI"
      MOTIF_STRING[ UInt8(SKIP_MOTIF) + 1 ] = "SE"
      MOTIF_STRING[ UInt8(ALTD_MOTIF) + 1 ] = "AD"
      MOTIF_STRING[ UInt8(ALTA_MOTIF) + 1 ] = "AA"
      MOTIF_STRING[ UInt8(NONE_MOTIF) + 1 ] = "na"

const MOTIF_TABLE = fill(NONE_MOTIF, 2^6 )
      # Alt TxStart  SL SL
      MOTIF_TABLE[ 0b000000 + 1 ] = TXST_MOTIF
      # Alt PolyA    RS RS
      MOTIF_TABLE[ 0b011011 + 1 ] = TXEN_MOTIF

      # Alt FirstEx  SL LS
      MOTIF_TABLE[ 0b000010 + 1 ] = ALTF_MOTIF
      # Alt LastEx   SR RS
      MOTIF_TABLE[ 0b001011 + 1 ] = ALTL_MOTIF

      # RetainedInt  LL RR
      MOTIF_TABLE[ 0b101110 + 1 ] = RETI_MOTIF

      # Skipped      RR LL
      MOTIF_TABLE[ 0b110101 + 1 ] = SKIP_MOTIF
      #              RR LR
      MOTIF_TABLE[ 0b110100 + 1 ] = SKIP_MOTIF
      #              LR LL
      MOTIF_TABLE[ 0b100101 + 1 ] = SKIP_MOTIF
      #              LR LR
      MOTIF_TABLE[ 0b100100 + 1 ] = SKIP_MOTIF

      # Alt Donor    LL LL
      MOTIF_TABLE[ 0b101101 + 1 ] = ALTD_MOTIF
      #              LL LR
      MOTIF_TABLE[ 0b101100 + 1 ] = ALTD_MOTIF

      # Alt Acceptor RR RR
      MOTIF_TABLE[ 0b110110 + 1 ] = ALTA_MOTIF
      #              LR RR
      MOTIF_TABLE[ 0b100110 + 1 ] = ALTA_MOTIF

Base.convert(::Type{EdgeMotif}, tup::Tuple{EdgeType,EdgeType}) = Base.convert(EdgeMotif, tup... )
Base.convert(::Type{EdgeMotif}, current::EdgeType, next::EdgeType) = MOTIF_TABLE[ (UInt8(current) << 3) | UInt8(next) + 1 ]

Base.convert{S <: AbstractString}(::Type{S}, edg::EdgeMotif ) = MOTIF_STRING[ UInt8(edg) + 1 ]

isobligate(  motif::EdgeMotif ) = motif != NONE_MOTIF && !( UInt8(motif) & 0b100 == 0b100 )
isaltsplice( motif::EdgeMotif ) = (UInt8(motif) & 0b110) == 0b110 

isspanning{I <: AbstractInterval, T <: Integer}( edge::I, node::T ) = edge.first < node < edge.last ? true : false
isconnecting{I <: AbstractInterval, T <: Integer}( edge::I, node::T ) = edge.first == node || edge.last == node ? true : false

# this is meant for short arrays when it is faster
# than using the overhead of a set
# DEPRECATED, IntSet is very efficient
function unique_push!{T}( arr::Vector{T}, el::T )
   if !( el in arr )
      push!( arr, el )
   end
end

type PsiPath
   psi::Float64
   count::Float64
   length::Float64
   nodes::IntSet
end

# This holds one or many sets of connected
# nodes + the count of the reads and the eff_len
# We also store the min and max node of all sets
# Note: It did occur to me that a cleaner solution
# is PsiGraph holding paths::Vector{PsiPath}
# but the current implementation was choosen for 
# faster memory access w.r.t. stride.
type PsiGraph
   psi::Vector{Float64}
   count::Vector{Float64}
   length::Vector{Float64}
   nodes::Vector{IntSet}
   min::NodeInt
   max::NodeInt
end

function hasintersect( a::IntSet, b::IntSet )
   seta,setb = length(a) > length(b) ? (b,a) : (a,b)
   for elem in seta
      if elem in setb
         return true
      end
   end
   false
end

function Base.in{T <: Integer}( i::T, pgraph::PsiGraph )
   for nodeset in pgraph.nodes
      (i in nodeset) && return true
   end
   return false
end

function reduce_graph!( pgraph::PsiGraph )
   i = 1
   @assert( length(pgraph.nodes) == length(pgraph.count) == length(pgraph.length) )
   while i < length(pgraph.nodes)
      if hasintersect( pgraph.nodes[i], pgraph.nodes[i+1] )
         for n in pgraph.nodes[i]
            push!( pgraph.nodes[i+1], n )
         end
         pgraph.count[i+1]  += pgraph.count[i]
         pgraph.length[i+1] += pgraph.length[i]
         splice!( pgraph.nodes, i )
         splice!( pgraph.count, i )
         splice!( pgraph.length, i )
         i -= 1
      end
      i += 1
   end
end

function Base.push!{I <: AbstractInterval}( pgraph::PsiGraph, edg::I; 
                                            value_bool=true, length=1.0 )
   push!( pgraph.count, value_bool ? edg.value : 0.0 )
   push!( pgraph.length, 1.0 )
   push!( pgraph.nodes, IntSet([edg.first, edg.last]) )
   reduce_graph!( pgraph )
   if edg.first < pgraph.min
      pgraph.min = edg.first
   end
   if edg.last > pgraph.max
      pgraph.max = edg.last
   end
end

function Base.push!{I <: AbstractInterval}( ppath::PsiPath, edg::I; 
                                            value_bool=true, length=1.0 )
   value_bool && (ppath.count += edg.value)
   ppath.length += length
   push!( ppath.nodes, edg.first )
   push!( ppath.nodes, edg.last  )
end

type AmbigCounts
   paths::Vector{NodeInt}
   prop::Vector{Float64}
   prop_sum::Float64
   multiplier::Float64
end

Base.(:(==))( a::AmbigCounts, b::AmbigCounts ) = a.paths == b.paths
Base.(:(==))( a::AmbigCounts, b::IntSet ) = ( a.paths == b )
Base.(:(==))( a::IntSet, b::AmbigCounts ) = ( a == b.paths )
Base.(:(==)){I <: Integer}( a::Vector{I}, b::IntSet ) = ( b == a )
function Base.(:(==)){I <: Integer}( iset::IntSet, ivec::Vector{I} )
   length(iset) != length(ivec) && return false
   for intv in ivec
      !(intv in iset) && return false
   end
   true
end

function Base.sum( vec::Vector{AmbigCounts} )
   sum = 0.0
   for i in 1:length(vec)
      sum += vec[i].multiplier
   end
   sum
end

function add_node_counts!( ambig::Vector{AmbigCounts}, ipath::PsiPath, 
                           egraph::PsiGraph, sgquant::SpliceGraphQuant )

   minv = min( egraph.min, first(ipath.nodes) )
   maxv = max( egraph.max, last(ipath.nodes)  )
   iset = IntSet()
   for n in minv:maxv # lets go through all possible nodes
      if n in ipath.nodes
         push!( iset, 1 )
         ipath.length += sgquant.leng[n]
      end
      for i in 1:length(egraph.nodes)
         if n in egraph.nodes[i]
            push!( iset, i+1 )
            egraph.length[i] += sgquant.leng[n]
         end
      end
      #=                                            =#
      length( iset ) == 0 && continue # unspliced node
      if length( iset ) == 1 # non-ambiguous node
         if first( iset ) == 1 # belongs to inclusion
            ipath.count += sgquant.node[n]
         else # belongs to exclusion
            idx = first( iset ) - 1
            egraph.count[idx] += sgquant.node[n]
         end
      else
         #check if there is already an entry for this set of paths
         # if so, just increment multiplier, if not, make new one
         exists = false 
         for am in ambig
            if iset == am
               am.multiplier += sgquant.node[n]
               exists = true
               break
            end
         end
         if !exists
            push!( ambig, AmbigCounts( collect(iset), 
                                       ones( length(iset) ) / length(iset), 
                                       1.0, sgquant.node[n] ) )
         end
      end
      empty!( iset ) # clean up
   end
end

function add_edge_counts!( ambig::Vector{AmbigCounts}, ipath::PsiPath,
                           egraph::PsiGraph, edges::Vector{IntervalValue} )
   iset = IntSet()
   for edg in edges
      if edg.first in ipath.nodes && edg.last in ipath.nodes
         push!(iset, 1)
         ipath.length += 1
      end
      for i in 1:length(egraph.nodes)
         if edg.first in egraph.nodes[i] && edg.last in egraph.nodes[i]
            push!( iset, i+1 )
            egraph.length[i] += 1
         end
      end
      #=                                             =#
      length( iset ) == 0 && continue
      if length( iset ) == 1
         if first( iset ) == 1 
            ipath.count += edg.value
         else
            idx = first( iset ) - 1
            egraph.count[idx] += edg.value
         end
      else
         exists = false
         for am in ambig
            if iset == am
               am.multiplier += edg.value
               exists = true
               break
            end
         end 
         if !exists
            push!( ambig, AmbigCounts( collect(iset),
                                       ones( length(iset) ) / length(iset),
                                       1.0, edg.value ) )
         end
      end
      empty!( iset )
   end
end

function _process_spliced( sg::SpliceGraph, sgquant::SpliceGraphQuant, 
                           node::NodeInt, motif::EdgeMotif )

   inc_path   = Nullable{PsiPath}()
   exc_graph  = Nullable{PsiGraph}()
   ambig_edge = Nullable{Vector{IntervalValue}}()
   ambig_cnt  = Nullable{Vector{AmbigCounts}}()

   for edg in intersect( sgquant.edge, (node, node) )
      if   isconnecting( edg, node )
         if isnull( inc_path )
            inc_path = Nullable(PsiPath( 0.0, 0.0, 0.0, IntSet() ))
         end   
         push!( inc_path.value, edg )
      elseif isspanning( edg, node )
         if isnull( exc_graph ) #don't allocate unless there is alt splicing
            exc_graph = Nullable(PsiGraph( Vector{Float64}(), Vector{Float64}(),
                                           Vector{Float64}(), Vector{IntSet}(),
                                           edg.first, edg.last ))
         end
         push!( exc_graph.value, edg )
      else
         error("Edge has to be connecting or spanning!!!!" )
      end
   end

   if !isnull( exc_graph ) && !isnull( inc_path )
      ambig_cnt = Nullable( Vector{AmbigCounts}() )
      # if the min or max of any exclusion set is different than the min/max
      # of the inclusion set we have a disjoint graph module and we can go
      # ahead and try to bridge nodes by extending with potentially ambiguous edges
      if !isnull( inc_path ) && ( get(exc_graph).min == first(get(inc_path).nodes) ||
                                  get(exc_graph).max == last(get(inc_path).nodes) )
         ambig_edge = extend_edges!( sgquant.edge, exc_graph.value, inc_path.value, ambig_edge, node )
         !isnull( ambig_edge ) && add_edge_counts!( ambig_cnt.value, inc_path.value, 
                                                    exc_graph.value, get(ambig_edge) )
      end
      add_node_counts!( ambig_cnt.value, inc_path.value, exc_graph.value, sgquant )
   end # end expanding module

  # lets finish up now.
   psi = Nullable{Float64}()
   if isnull( exc_graph ) # no spanning edge
      # check if we have both inclusion edges represented, or one if alt 5'/3'
      if !isnull( inc_path ) && ((get(inc_path).length >= 2 && get(inc_path).count >= 1) || 
                                 (get(inc_path).length >= 1 && get(inc_path).count >= 1  && 
                                  isaltsplice(motif)))
          psi = Nullable( 0.99 ) #&& likelihood_ci( psi, inc_cnt, z=1.64 )
      else
         # NA is ignored.
      end
   else # there is skipping
      if isnull( inc_path )
         if sum( get(exc_graph).count ) >= 1
             psi = Nullable( 0.0 ) #&& likelihood_ci( psi, exc_graph.count, z=1.64 )
         else
            # NA
         end
      else
         get(exc_graph).psi = zeros( length( get(exc_graph).count ) )
         calculate_psi!( inc_path.value, exc_graph.value, [get(inc_path).count; get(exc_graph).count] )
         println( "$inc_path\n$exc_graph\n$ambig_cnt\n\n" )
         it = rec_spliced_em!( inc_path.value, exc_graph.value, ambig_cnt.value, sig=4 )
         println( "$inc_path\n$exc_graph\n$ambig_cnt\n\n" )
         psi = Nullable( get(inc_path).psi )
      end
   end
   psi,inc_path,exc_graph,ambig_cnt
end

function extend_edges!{K,V}( edges::IntervalMap{K,V}, pgraph::PsiGraph, ipath::PsiPath,
                             ambig_edge::Nullable{Vector{IntervalValue}}, node::NodeInt )

   minv = min( pgraph.min, first(ipath.nodes) )
   maxv = max( pgraph.max, last(ipath.nodes)  )
   idx = node - 1
   shouldpush = false
   while idx > minv # <----- node iter left
      #iterate through local edges to the left
      for edg in intersect( edges, (idx,idx) )
         # if this is a new edge and connects to our idx from the left
         if isconnecting( edg, idx ) && edg.last == idx && edg.first >= minv
            shouldpush = false
            if edg.last in pgraph
               push!( pgraph, edg, value_bool=false )
               shouldpush = true
            end
            if edg.last in ipath.nodes
               push!( ipath.nodes, edg.first )
               shouldpush = true
            end
            if shouldpush
               if isnull( ambig_edge )
                  ambig_edge = Nullable( Vector{IntervalValue}() )
               end
               push!( ambig_edge.value, edg )
            end
         end
      end
      idx -= 1
   end

   idx = node + 1
   while idx < maxv # iter right node ------->
      # iterate through local edges to the right
      for edg in intersect( edges, (idx,idx) )
         # if this is a new edge and connects to our idx from the left
         if isconnecting( edg, idx ) && edg.first == idx && edg.last <= maxv
            shouldpush = false
            if edg.last in pgraph
               push!( pgraph, edg, value_bool=false )
               shouldpush = true
            end
            if edg.last in ipath.nodes
               push!( ipath.nodes, edg.first )
               shouldpush = true
            end
            if shouldpush
               if isnull( ambig_edge )
                  ambig_edge = Nullable( Vector{IntervalValue}() )
               end
               push!( ambig_edge.value, edg )
            end
         end
      end
      idx += 1
   end

   ambig_edge
end

function process_events( outfile, lib::GraphLib, anno::Refset, graphq::GraphLibQuant )
   io = open( outfile, "w" )
   stream = ZlibDeflateOutputStream(io)
   for g in 1:length(lib.graphs)
      name = lib.names[g]
      chr,strand = anno.geneset[ name ].info
      println(STDERR, "$g, $name, $chr, $strand" )
      _process_events( stream, lib.graphs[g], graphq.quant[g], (name,chr,strand) )
   end
   close(stream)
   close(io)
end

typealias GeneMeta Tuple{Genename, Seqname, Char}
typealias BufOut BufferedStreams.BufferedOutputStream

function _process_events( io::BufOut, sg::SpliceGraph, sgquant::SpliceGraphQuant, info::GeneMeta )
   # Heres the plan:
   # step through sets of edges, look for edge motifs, some are obligate calculations
   # others we only calculate psi if there is an alternative edge
   # Then calculate bias of the event 
   if length( sg.edgetype ) <= 2
      return 0
   end
   for i in 1:length(sg.edgetype)-1
      motif = convert(EdgeMotif, sg.edgetype[i], sg.edgetype[i+1] )
      motif == NONE_MOTIF && continue
      if isobligate( motif ) # is utr event
          
      else  # is a spliced node
         psi,inc,exc,ambig = _process_spliced( sg, sgquant, convert(NodeInt, i) , motif )
         if !isnull( psi )
            ambig_cnt = isnull( ambig ) ? 0.0 : sum( ambig.value )
            output_psi( io, get(psi), inc, exc, ambig_cnt, motif, sg, i, info )
         end
      end  
   end
end

tab_write{S <: AbstractString}( io::BufOut, str::S ) = (write( io, str ); write( io, '\t'  ))
function coord_write( io::BufOut, chr, first, last )
   write( io, chr   )
   write( io, ':'   )
   write( io, string(first) )
   write( io, '-'   )
   write( io, string(last)  )
end
function coord_write( io::BufOut, chr, first, last, strand )
   coord_write( io, chr, first, last )
   write( io, ':'    )
   write( io, strand )
end

function output_psi( io::BufOut, psi::Float64, inc::Nullable{PsiPath}, exc::Nullable{PsiGraph}, 
                     ambig::Float64, motif::EdgeMotif, sg::SpliceGraph, node::Int, 
                     info::GeneMeta )
    
   # gene
   tab_write( io, info[1] )
   # coordinate
   coord_write( io, info[2], sg.nodecoord[node], sg.nodecoord[node]+sg.nodelen[node]-1, info[3] )
   write( io, '\t' )
   # event_type
   tab_write( io, convert(ASCIIString, motif) )
   # psi
   tab_write( io, string(psi) )
   
   if !isnull( inc ) && !isnull( exc )
      tab_write( io, string(get(inc).nodes) )

      for i in 1:length(get(exc).nodes)-1
         write( io, string( get(exc).nodes[i] ) )
         write( io, "," )
      end
      write( io, string( get(exc).nodes[end] ) )
   end

   write( io, '\n' )
end

function Base.unsafe_copy!{T <: Number}( dest::Vector{T}, src::Vector{T}; indx_shift=0 )
   for i in 1:length(src)
      dest[i+indx_shift] = src[i]
   end
end

function calculate_psi!( ipath::PsiPath, egraph::PsiGraph, counts::Vector{Float64}; sig=0 )
   ipath.psi = counts[1] / ipath.length
   for i in 1:length(egraph.psi)
      egraph.psi[i] = counts[i+1] / egraph.length[i]
   end
   cnt_sum = sum( egraph.psi ) + ipath.psi
   if sig > 0
      @fastmath ipath.psi = signif( ipath.psi / cnt_sum, sig )
      for i in 1:length( egraph.psi )
         @fastmath egraph.psi[i] = signif( egraph.psi[i] / cnt_sum, sig )
      end 
   else
      @fastmath ipath.psi /= cnt_sum
      for i in 1:length( egraph.psi )
         @fastmath egraph.psi[i] /= cnt_sum
      end
   end
end

function rec_spliced_em!( ipath::PsiPath, egraph::PsiGraph, 
                          ambig::Vector{AmbigCounts};
                          inc_temp::Float64=0.0,
                          exc_temp::Vector{Float64}=zeros(1+length(egraph.count)),
                          count_temp::Vector{Float64}=ones(1+length(egraph.count)),
                          it=1, max=1000, sig=0 )

   count_temp[1] = ipath.count
   unsafe_copy!( count_temp, egraph.count, indx_shift=1 )
   inc_temp = ipath.psi
   unsafe_copy!( exc_temp, egraph.psi )
   
   for ac in ambig
      idx = 1
      if it > 1 # maximization
         ac.prop_sum = 0.0
         for p in ac.paths
            prev_psi = p == 1 ? ipath.psi : egraph.psi[p-1]
            ac.prop[idx] = prev_psi
            ac.prop_sum += prev_psi
            idx += 1
         end
      end

      idx = 1    
      for p in ac.paths
         @fastmath prop = ac.prop[idx] / ac.prop_sum
         ac.prop[idx] = prop
         @fastmath count_temp[p] += prop * ac.multiplier
         idx += 1
      end
   end

   calculate_psi!( ipath, egraph, count_temp, sig=sig ) # expectation

   if (inc_temp != ipath.psi || egraph.psi != exc_temp) && it < max
      it = rec_spliced_em!( ipath, egraph, ambig, 
                            inc_temp=inc_temp, exc_temp=exc_temp, count_temp=count_temp,
                            it=it+1, max=max, sig=sig )
   end
   it
end
