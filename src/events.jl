
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

isobligate(  motif::EdgeMotif ) = motif != NONE_MOTIF && !( UInt8(motif) & 0b100 == 0b100 )
isaltsplice( motif::EdgeMotif ) = (UInt8(motif) & 0b110) == 0b110 

function process_sgquant( lib::GraphLib, graphq::GraphLibQuant )
   for g in 1:length(lib.graphs)
      process_events!( lib.graphs[g], graphq.quant[g] )
   end
end

isspanning{I <: AbstractInterval}( edge::I, node::Coordint ) = edge.first < node < edge.last ? true : false
isconnecting{I <: AbstractInterval}( edge::I, node::Coordint ) = edge.first == node || edge.last == node ? true : false

# this is meant for short arrays when it is faster
# than using the overhead of a set
# DEPRECATED, IntSet is very efficient
function unique_push!{T}( arr::Vector{T}, el::T )
   if !( el in arr )
      push!( arr, el )
   end
end

# This holds one or many sets of connected
# nodes + the count of the reads and the eff_len
# We also store the min and max node of all sets
type PsiGraph
   count::Vector{Float64}
   length::Vector{Float64}
   nodes::Vector{IntSet}
   min::Coordint
   max::Coordint
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

function Base.in( i::Int, pgraph::PsiGraph )
   for nodeset in pgraph.nodes
      (i in nodeset) && return true
   end
   return false
end

function reduce_graph!( pgraph::PsiGraph )
   i = 1
   while i < length(pgraph.nodes)
      if hasintersect( pgraph.nodes[i], pgraph.nodes[i+1] )
         pgraph.nodes[i+1]   = union( pgraph.nodes[i], pgraph.nodes[i+1] )
         pgraph.count[i+1]  += pgraph.count[i]
         pgraph.length[i+1] += pgraph.length[i]
         shift!( pgraph.nodes  )
         shift!( pgraph.count  )
         shift!( pgraph.length )
         i -= 1
      end
      i += 1
   end
end

function Base.push!{I <: AbstractInterval}( pgraph::PsiGraph, edg::I; value_bool=true, length=1.0 )
   value_bool && push!( pgraph.count, edg.value )
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

function _process_spliced_pg( sg::SpliceGraph, sgquant::SpliceGraphQuant, node::Coordint, motif::EdgeMotif, eff_len::Int )
   inc_cnt = 0.0
   inc_len = 0.0
   inc_set = IntSet()
   exc_graph  = Nullable{PsiGraph}()
   ambig_edge = Nullable{Vector{IntervalValue}}()

   for edg in intersect( sgquant.edge, (node, node) )
      if   isconnecting( edg, node )
         inc_cnt += edg.value
         inc_len += 1
         push!( inc_set, edg.first )
         push!( inc_set, edg.last  )
      elseif isspanning( edg, node )
         if isnull( exc_graph ) #don't allocate unless there is alt splicing
            exc_graph = Nullable(PsiGraph( Vector{Float64}(), 
                                           Vector{Float64}(), 
                                           Vector{IntSet}(), edg.first, edg.last ))
         end
         push!( get(exc_graph), edg )
      else
         error("Edge has to be connecting or spanning!!!!" )
      end
   end

   # if the min or max of any exclusion set is different than the min/max
   # of the inclusion set we have a disjoint graph module and we can go
   # ahead and try to bridge nodes by iteratively adding ambiguous edges
   if !isnull( exc_graph ) && ( get(exc_graph).min == first(inc_set) ||
                                get(exc_graph).max == last(inc_set) )

      ambig_edge,inc_len = extend_edges!( sgquant.edge, get(exc_graph),
                                          inc_set, inc_len, ambig_edge )
   end # end expanding module

   # now we need to make ambiguous counts, and then do EM

end

function extend_edges!( edges::IntervalMap, pgraph::PsiGraph, 
                        inc_set::IntSet, inc_len::Float64,
                        ambig_edge::Nullable{Vector{IntervalValue}} )

   min = min( pgraph.min, first(inc_set) )
   max = max( pgraph.max, last(inc_set)  )
   idx = node - 1
   shouldpush = false
   while idx > min # <----- node iter left
      #iterate through local edges to the left
      for edg in intersect( edges, (idx,idx) )
         # if this is a new edge and connects to our idx from the left
         if isconnecting( edg, idx ) && edg.last == idx && edg.first >= min
            shouldpush = false
            if edg.last in pgraph
               push!( pgraph, edg, value_bool=false )
               shouldpush = true
            end
            if edg.last in inc_set
               push!( inc_set, edg.first )
               inc_len += 1
               shouldpush = true
            end
            if shouldpush
               if isnull( ambig_edge )
                  ambig_edge = Nullable( Vector{IntervalValue}() )
               end
               push!( ambig_edge, edg )
            end
         end
      end
      idx -= 1
   end

   idx = node + 1
   while idx < max # iter right node ------->
      # iterate through local edges to the right
      for edg in intersect( edges, (idx,idx) )
         # if this is a new edge and connects to our idx from the left
         if isconnecting( edg, idx ) && edg.first == idx && edg.last <= min
            shouldpush = false
            if edg.last in pgraph
               push!( pgraph, edg, value_bool=false )
               shouldpush = true
            end
            if edg.last in inc_set
               push!( inc_set, edg.first )
               inc_len += 1
               shouldpush = true
            end
            if shouldpush
               if isnull( ambig_edge )
                  ambig_edge = Nullable( Vector{IntervalValue}() )
               end
               push!( ambig_edge, edg )
            end
         end
      end
      idx += 1
   end

   ambig_edge, inc_len
end


function _process_spliced( sg::SpliceGraph, sgquant::SpliceGraphQuant, node::Coordint, motif::EdgeMotif, eff_len::Int )
   inc_cnt = 0.0
   inc_len = 0.0
   inc_set = IntSet()
   exc_cnt = Nullable{Vector{Float64}}()
   exc_len = Nullable{Vector{Float64}}()
   exc_set = Nullable{Vector{IntSet}}()

   # push initial graph structure for inclusion/exclusion-set
   for edg in intersect( sgquant.edge, (node, node) )
      if   isconnecting( edg, node )
         inc_cnt += edg.value
         inc_len += 1
         push!( inc_set, edg.first )
         push!( inc_set, edg.last  )
      elseif isspanning( edg, node )
         if isnull( exc_set ) #don't allocate unless there is alt splicing
            exc_cnt = Nullable(Vector{Float64}())
            exc_len = Nullable(Vector{Float64}())
            exc_set = Nullable(Vector{IntSet}())
         end
         push!( get(exc_cnt), edg.value )
         push!( get(exc_len), 1.0 )
         push!( get(exc_set), IntSet([edg.first, edg.last]) )
         reduce_graph!( get(exc_set), get(exc_cnt), get(exc_len) )
      else
         error("Edge has to be connecting or spanning!!!!" )
      end
   end

   # if the min or max of any exclusion set is greater than the min/max
   # of the inclusion set we have a disjoint graph module and we can go
   # ahead and try to bridge nodes by iteratively adding ambiguous edges
   if !isnull( exc_set )
      while 
   end 

   # lets finish up now.
   if exc_len == 0 # no spanning edge
      # check if we have both inclusion edges represented, or one if alt 5'/3'
      if (inc_len >= 2 && inc_cnt >= 1) || (inc_len >= 1 && inc_cnt >= 1  && isaltsplice(motif)) 
         # psi = 0.99 && likelihood_ci( psi, inc_cnt, z=1.64 )
      else
         # NA is ignored.
      end
   else # do EM
      psi = rec_spliced_em!(  )
   end
   #conf_int = likelihood_ci()
end


function output_psi{F <: AbstractFloat}( icnt::F, ecnt::F, ilen::F, elen::F, 
                                         nodes::Vector{Coordint}, node::Coordint )
   
end

function process_events( sg::SpliceGraph, sgquant::SpliceGraphQuant, eff_len::Int )
   # Heres the plan:
   # step through sets of edges, look for edge motifs, some are obligate calculations
   # others we only calculate psi if there is an alternative edge
   # Then calculate bias of the event 
   if length( sg.edges ) <= 2
      return 0
   end
   for i in 1:length(sg.edges)-1
      motif = convert(EdgeMotif, sg.edges[i], sg.edges[i+1] )
      motif == NONE_MOTIF && continue
      if isobligate( motif ) # is utr event
          
      else  # is a spliced node
         _process_spliced( sg, sgquant, sg.edges[i], motif, eff_len )
      end  
   end
end


