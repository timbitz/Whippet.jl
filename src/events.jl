
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
      MOTIF_STRING[ UInt8(SKIP_MOTIF) + 1 ] = "CE"
      MOTIF_STRING[ UInt8(ALTD_MOTIF) + 1 ] = "AD"
      MOTIF_STRING[ UInt8(ALTA_MOTIF) + 1 ] = "AA"
      MOTIF_STRING[ UInt8(NONE_MOTIF) + 1 ] = "na"

const MOTIF_TABLE = fill(NONE_MOTIF, 2^6 )
      # Alt TxStart  SL SL
      MOTIF_TABLE[ 0b000000 + 1 ] = TXST_MOTIF
      # Alt TxStart  LS SL
      MOTIF_TABLE[ 0b010000 + 1 ] = TXST_MOTIF

      # Alt PolyA    RS RS
      MOTIF_TABLE[ 0b011011 + 1 ] = TXEN_MOTIF
      # Alt PolyA    RS SR
      MOTIF_TABLE[ 0b011001 + 1 ] = TXEN_MOTIF

      # Alt FirstEx  LS LL
      MOTIF_TABLE[ 0b010101 + 1 ] = ALTF_MOTIF
      # Alt FirstEx  LS LR
      MOTIF_TABLE[ 0b010100 + 1 ] = ALTF_MOTIF
      # Alt FirstEx  LS LS
      MOTIF_TABLE[ 0b010010 + 1 ] = ALTF_MOTIF

      # Alt LastEx   RR SR
      MOTIF_TABLE[ 0b110001 + 1 ] = ALTL_MOTIF
      # Alt LastEx   LR SR
      MOTIF_TABLE[ 0b100001 + 1 ] = ALTL_MOTIF
      # Alt LastEx   SR SR 
      MOTIF_TABLE[ 0b001001 + 1 ] = ALTL_MOTIF

      # RetainedInt  LL RR
      MOTIF_TABLE[ 0b101110 + 1 ] = RETI_MOTIF

      # Skipped      RR LL
      MOTIF_TABLE[ 0b110101 + 1 ] = SKIP_MOTIF
      #              RR LR
      MOTIF_TABLE[ 0b110100 + 1 ] = SKIP_MOTIF
      #              RR LS
      MOTIF_TABLE[ 0b110010 + 1 ] = SKIP_MOTIF
     
      #              LR LL
      MOTIF_TABLE[ 0b100101 + 1 ] = SKIP_MOTIF
      #              LR LR
      MOTIF_TABLE[ 0b100100 + 1 ] = SKIP_MOTIF
      #              LR LS
      MOTIF_TABLE[ 0b100010 + 1 ] = SKIP_MOTIF

      #              SR LL
      MOTIF_TABLE[ 0b001101 + 1 ] = SKIP_MOTIF
      #              SR LR
      MOTIF_TABLE[ 0b001100 + 1 ] = SKIP_MOTIF
      #              SR LS
      MOTIF_TABLE[ 0b001010 + 1 ] = SKIP_MOTIF


      # Alt Donor    LL LL
      MOTIF_TABLE[ 0b101101 + 1 ] = ALTD_MOTIF
      #              LL LR
      MOTIF_TABLE[ 0b101100 + 1 ] = ALTD_MOTIF
      #              LL LS
      MOTIF_TABLE[ 0b101010 + 1 ] = ALTD_MOTIF

      # Alt Acceptor RR RR
      MOTIF_TABLE[ 0b110110 + 1 ] = ALTA_MOTIF
      #              LR RR
      MOTIF_TABLE[ 0b100110 + 1 ] = ALTA_MOTIF
      #              SR RR
      MOTIF_TABLE[ 0b001110 + 1 ] = ALTA_MOTIF

Base.convert(::Type{EdgeMotif}, tup::Tuple{EdgeType,EdgeType}) = Base.convert(EdgeMotif, tup... )
Base.convert(::Type{EdgeMotif}, current::EdgeType, next::EdgeType) = MOTIF_TABLE[ (UInt8(current) << 3) | UInt8(next) + 1 ]

Base.convert{S <: AbstractString}(::Type{S}, edg::EdgeMotif ) = MOTIF_STRING[ UInt8(edg) + 1 ]

#isobligate(  motif::EdgeMotif ) = motif != NONE_MOTIF && !( UInt8(motif) & 0b100 == 0b100 )
isobligate( motif::EdgeMotif ) = (0 <= UInt8(motif) <= 1)
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

# Deprecated 
type PsiPath
   psi::Float64
   count::Float64
   length::Float64
   nodes::IntSet
end # TODO: Delete PsiPath.


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

Base.sum( pgraph::PsiGraph ) = sum( pgraph.count )
Base.sum( ngraph::Nullable{PsiGraph} ) = isnull( ngraph ) ? 
                                         0.0 : convert(Float64, sum( ngraph.value ))

function hasintersect( a::IntSet, b::IntSet )
   seta,setb = length(a) > length(b) ? (b,a) : (a,b)
   for elem in seta
      if elem in setb
         return true
      end
   end
   false
end

function hasintersect_terminal( a::IntSet, b::IntSet )
   if first(a) == last(b) || first(b) == last(a)
      return true
   else
      return false
   end
end

function Base.in{T <: Integer}( i::T, pgraph::PsiGraph )
   for nodeset in pgraph.nodes
      (i in nodeset) && return true
   end
   return false
end

# This looks for an edge in an IntSet... The two nodes in the edge
# must be adjacent in the IntSet to be true.
function Base.in{I <: IntervalValue}( edge::I, iset::IntSet )
   s = start(iset)
   while !done( iset, s )
      v1,s = next( iset, s )
      if v1 == edge.first
         v2,s = next( iset, s )
         if v2 == edge.last
            return true
         end
      end
   end
   false
end

complexity( one::PsiGraph, two::PsiGraph ) = complexity(length(one.nodes) + length(two.nodes))
complexity( one::PsiGraph ) = complexity(length(one.nodes))
# This function calculates complexity of a splicing event
# as the ceil log2 number of paths through the graph
complexity( num_paths::Int ) = @fastmath Int(ceil(log2( num_paths )))

# Build minimal set of paths to explain graph
# Explanation: If we have 3 edges for example 1-2, 2-3, and 1-3
# We have to create 2 paths, 1-2-3, and 1-3, since 1-3 can't
# possibly contain the second node.
# This function will do this with any number of edges or complex
# paths through a graph and is the essence of how Whippet quantification
# works.. these graph paths are then assigned unambiguous counts, while
# ambiguous counts are then assigned to each graph path using the EM algorithm
function reduce_graph( pgraph::PsiGraph )
   newgraph = PsiGraph( Vector{Float64}(), Vector{Float64}(),
                        Vector{Float64}(), Vector{IntSet}(),
                        pgraph.min, pgraph.max )
   used = IntSet()
   for i in 1:length(pgraph.nodes)
      if i < length(pgraph.nodes)
         for j in i:length(pgraph.nodes)
            if hasintersect_terminal( pgraph.nodes[i], pgraph.nodes[j] )
               push!( newgraph.nodes,  union( pgraph.nodes[i], pgraph.nodes[j] ) )
               push!( newgraph.length, pgraph.length[i] + pgraph.length[j] )
               push!( newgraph.count,  pgraph.count[i] + pgraph.count[j] )
               push!( used, i )
               push!( used, j )
            end
         end
      end
      if !( i in used )
         push!( newgraph.nodes, pgraph.nodes[i] )
         push!( newgraph.length, pgraph.length[i] )
         push!( newgraph.count, pgraph.count[i] )
         push!( used, i )
      end
   end
   ischanging = true
   while ischanging
      ischanging = false
      for i in 1:length(pgraph.nodes)
         for j in 1:length(newgraph.nodes)
            if hasintersect_terminal( pgraph.nodes[i], newgraph.nodes[j] )
               for n in pgraph.nodes[i]
                  push!( newgraph.nodes[j], n )
               end
               newgraph.count[j] += pgraph.count[i]
               newgraph.length[j] += pgraph.length[i]
               ischanging = true 
            end
         end
      end
   end
   i = 1
   if length(newgraph.nodes) > 1
     while i < length(newgraph.nodes)
        if newgraph.nodes[i] in newgraph.nodes[i+1:end]
           splice!( newgraph.nodes,  i )
           splice!( newgraph.length, i )
           splice!( newgraph.count,  i )
           i -= 1
        end
        i += 1
     end
   end
   newgraph
end

function Base.push!{I <: AbstractInterval}( pgraph::PsiGraph, edg::I; 
                                            value_bool=true, length=1.0 )
   push!( pgraph.count, value_bool ? edg.value : 0.0 )
   push!( pgraph.length, value_bool ? length : 0.0 )
   push!( pgraph.nodes, IntSet([edg.first, edg.last]) )
   #reduce_graph_terminal!( pgraph )
   if edg.first < pgraph.min
      pgraph.min = edg.first
   end
   if edg.last > pgraph.max
      pgraph.max = edg.last
   end
end

function Base.push!{I <: AbstractInterval}( ppath::PsiPath, edg::I; 
                                            value_bool=true, length=1.0 )
   ppath.count  += (value_bool ? edg.value : 0.0)
   ppath.length += (value_bool ? length : 0.0)
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

Base.sum( nvec::Nullable{Vector{AmbigCounts}} ) = isnull( nvec ) ? 0.0 :
                                                  sum( nvec.value )
function Base.sum( vec::Vector{AmbigCounts} )
   sum = 0.0
   for i in 1:length(vec)
      sum += vec[i].multiplier
   end
   sum
end

# add_node_counts! for single pgraph (used in tandem_utrs)
function add_node_counts!( ambig::Vector{AmbigCounts}, pgraph::PsiGraph,
                           sgquant::SpliceGraphQuant, bias::Float64 )

   iset = IntSet()
   for n in pgraph.min:pgraph.max # lets go through all possible nodes
      for i in 1:length(pgraph.nodes)
         if n in pgraph.nodes[i]
            push!( iset, i )
            pgraph.length[i] += 1 # sgquant.leng[n]
         end
      end

      #=                                            =#
      length( iset ) == 0 && continue # unspliced node
      if length( iset ) == 1 # non-ambiguous node
         if first( iset ) <= length(pgraph.nodes) # belongs to inclusion
            idx = first( iset )
            @fastmath pgraph.count[idx] += sgquant.node[n] * bias / sgquant.leng[n]
         end
      else
         #check if there is already an entry for this set of paths
         # if so, just increment multiplier, if not, make new one
         exists = false
         for am in ambig
            if iset == am
               @fastmath am.multiplier += sgquant.node[n] * bias / sgquant.leng[n]
               exists = true
               break
            end
         end
         if !exists && sgquant.node[n] > 0
            push!( ambig, AmbigCounts( collect(iset),
                                       ones( length(iset) ) / length(iset),
                                       1.0, sgquant.node[n] * bias / sgquant.leng[n] ) )
         end
      end
      empty!( iset ) # clean up
   end
end

# add_node_counts! for two graphs, inc_graph and exc_graph (used in spliced events)
function add_node_counts!( ambig::Vector{AmbigCounts}, igraph::PsiGraph, 
                           egraph::PsiGraph, sgquant::SpliceGraphQuant, 
                           bias::Float64 )

   minv = min( egraph.min, igraph.min )
   maxv = max( egraph.max, igraph.max )
   iset = IntSet()
   for n in minv:maxv # lets go through all possible nodes
      for i in 1:length(igraph.nodes)
         if n in igraph.nodes[i]
            push!( iset, i )
            igraph.length[i] += 1 # sgquant.leng[n]
         end
      end
      for i in 1:length(egraph.nodes)
         if n in egraph.nodes[i]
            push!( iset, i+length(igraph.nodes) )
            egraph.length[i] += 1 # sgquant.leng[n]
         end
      end
      #=                                            =#
      length( iset ) == 0 && continue # unspliced node
      if length( iset ) == 1 # non-ambiguous node
         if first( iset ) <= length(igraph.nodes) # belongs to inclusion
            idx = first( iset )
            @fastmath igraph.count[idx] += sgquant.node[n] * bias / sgquant.leng[n] 
         else # belongs to exclusion
            idx = first( iset ) - length(igraph.nodes)
            @fastmath egraph.count[idx] += sgquant.node[n] * bias / sgquant.leng[n]
         end
      else
         #check if there is already an entry for this set of paths
         # if so, just increment multiplier, if not, make new one
         exists = false 
         for am in ambig
            if iset == am
               @fastmath am.multiplier += sgquant.node[n] * bias / sgquant.leng[n]
               exists = true
               break
            end
         end
         if !exists && sgquant.node[n] > 0
            push!( ambig, AmbigCounts( collect(iset), 
                                       ones( length(iset) ) / length(iset), 
                                       1.0, sgquant.node[n] * bias / sgquant.leng[n] ) )
         end
      end
      empty!( iset ) # clean up
   end
end

# add_edge_counts! for one psigraph, (used in tandem utr)
function add_edge_counts!( ambig::Vector{AmbigCounts}, pgraph::PsiGraph,
                           sgquant::SpliceGraphQuant )
   iset = IntSet()
   for edg in sgquant.edge
      for i in 1:length(pgraph.nodes)
         if edg in pgraph.nodes[i]
            push!(iset, i)
            pgraph.length[i] += 1
         end
      end
      #=                                             =#
      length( iset ) == 0 && continue
      if length( iset ) == 1
         idx = first( iset )
         pgraph.count[idx] += edg.value
      else
         exists = false
         for am in ambig
            if iset == am
               am.multiplier += edg.value
               exists = true
               break
            end
         end
         if !exists && edg.value > 0
            push!( ambig, AmbigCounts( collect(iset),
                                       ones( length(iset) ) / length(iset),
                                       1.0, edg.value ) )
         end
      end
      empty!( iset )
   end
end

# add_edge_counts! for two graphs, inc_graph and exc_graph! (used in spliced events)
function add_edge_counts!( ambig::Vector{AmbigCounts}, igraph::PsiGraph,
                           egraph::PsiGraph, edges::Vector{IntervalValue} )
   iset = IntSet()
   for edg in edges
      for i in 1:length(igraph.nodes)
         if edg in igraph.nodes[i]
            push!(iset, i)
            igraph.length[i] += 1
         end
      end
      for i in 1:length(egraph.nodes)
         if edg in egraph.nodes[i]
            push!( iset, i+length(igraph.nodes) )
            egraph.length[i] += 1
         end
      end
      #=                                             =#
      length( iset ) == 0 && continue
      if length( iset ) == 1
         if first( iset ) <= length(igraph.nodes)
            idx = first( iset )
            igraph.count[idx] += edg.value
         else
            idx = first( iset ) - length(igraph.nodes)
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
         if !exists && edg.value > 0
            push!( ambig, AmbigCounts( collect(iset),
                                       ones( length(iset) ) / length(iset),
                                       1.0, edg.value ) )
         end
      end
      empty!( iset )
   end
end

function build_utr_graph( nodes::IntSet, motif::EdgeMotif, sgquant::SpliceGraphQuant )

   utr_graph = Nullable(PsiGraph( Vector{Float64}(), Vector{Float64}(),
                                  Vector{Float64}(), Vector{IntSet}(),
                                  first(nodes), last(nodes) ))
   curset = IntSet()

   while length(nodes) > 0
      if motif == TXST_MOTIF
         curval = pop!(nodes)
      else
         curval = shift!(nodes)
      end
      push!( curset, curval )
      push!( get(utr_graph).nodes, copy(curset) )
      push!( get(utr_graph).psi, 0.0 )
      push!( get(utr_graph).length, 0.0 )
      push!( get(utr_graph).count, 0.0 )
   end

   utr_graph
end

function _process_tandem_utr( sg::SpliceGraph, sgquant::SpliceGraphQuant,
                              node::NodeInt, motif::EdgeMotif )
   
   utr_graph  = Nullable{PsiGraph}()
   ambig_cnt  = Nullable{Vector{AmbigCounts}}()

   used_node  = IntSet()
   total_cnt  = 0.0

   if motif == TXEN_MOTIF
      # add node directly upstream (CE next to tandemUTR)
      push!( used_node, node-1 )
      total_cnt += sgquant.node[node-1]
   end

   i = node
   curmotif = motif
   while curmotif == motif
      push!( used_node, i )
      total_cnt += sgquant.node[i]
      if i > 1
         interv = Interval{Exonmax}( i-1, i )
         total_cnt += get( sgquant.edge, interv, IntervalValue(0,0,0.0) ).value
      end
      i += 1
      if i < length(sg.edgetype)
         curmotif = convert(EdgeMotif, sg.edgetype[i], sg.edgetype[i+1] )
      else
         break
      end
   end

   if motif == TXST_MOTIF
      # add node directly downstream (CE next to tandemUTR)
      push!( used_node, i )
      total_cnt += sgquant.node[i]
      interv = Interval{Exonmax}( i-1, i )
      total_cnt += get( sgquant.edge, interv, IntervalValue(0,0,0.0) ).value
   end

   psi = Nullable{Float64}()
   if total_cnt > 2
      ambig_cnt = Nullable( Vector{AmbigCounts}() )
      utr_graph = build_utr_graph( used_node, motif, sgquant )
      add_node_counts!( ambig_cnt.value, utr_graph.value, sgquant, 1.0 )
      add_edge_counts!( ambig_cnt.value, utr_graph.value, sgquant )
      it = rec_tandem_em!( utr_graph.value, ambig_cnt.value, sig=4 )
      psi = Nullable( get(utr_graph).psi )
   end

   psi,utr_graph,ambig_cnt
end

function _process_spliced( sg::SpliceGraph, sgquant::SpliceGraphQuant, 
                           node::NodeInt, motif::EdgeMotif, bias::Float64, isnodeok::Bool )

   inc_graph  = Nullable{PsiGraph}()
   exc_graph  = Nullable{PsiGraph}()
   ambig_edge = Nullable{Vector{IntervalValue}}()
   ambig_cnt  = Nullable{Vector{AmbigCounts}}()

   connecting_val = 0.0 
   spanning_val   = 0.0

   for edg in intersect( sgquant.edge, (node, node) )
      if   isconnecting( edg, node )
         if isnull( inc_graph )
            inc_graph = Nullable(PsiGraph( Vector{Float64}(), Vector{Float64}(),
                                           Vector{Float64}(), Vector{IntSet}(),
                                           edg.first, edg.last ))
         end   
         if isnull( ambig_edge )
            ambig_edge = Nullable( Vector{IntervalValue}() )
         end
         push!( inc_graph.value, edg, value_bool=false )
         push!( ambig_edge.value, edg )
         connecting_val += edg.value
      elseif isspanning( edg, node )
         if isnull( exc_graph ) #don't allocate unless there is alt splicing
            exc_graph = Nullable(PsiGraph( Vector{Float64}(), Vector{Float64}(),
                                           Vector{Float64}(), Vector{IntSet}(),
                                           edg.first, edg.last ))
         end
         if isnull( ambig_edge )
            ambig_edge = Nullable( Vector{IntervalValue}() )
         end
         push!( exc_graph.value, edg, value_bool=false )
         push!( ambig_edge.value, edg )
         spanning_val += edg.value
      else
         error("Edge has to be connecting or spanning!!!!" )
      end
   end

   if !isnull( exc_graph ) && !isnull( inc_graph )
      ambig_cnt = Nullable( Vector{AmbigCounts}() )

      # try to bridge nodes by extending with potentially ambiguous edges
      ambig_edge = extend_edges!( sgquant.edge, exc_graph.value, inc_graph.value, ambig_edge, node )

      inc_graph = Nullable( reduce_graph( inc_graph.value ) )
      exc_graph = Nullable( reduce_graph( exc_graph.value ) )

      add_edge_counts!( ambig_cnt.value, inc_graph.value, 
                        exc_graph.value, get(ambig_edge) )

      if isnodeok
         add_node_counts!( ambig_cnt.value, inc_graph.value, exc_graph.value, sgquant, bias )
      end
   end # end expanding module

  # lets finish up now.
   psi = Nullable{Float64}()
   if isnull( exc_graph ) # no spanning edge
      # check if we have both inclusion edges represented, or one if alt 5'/3'
      if !isnull( inc_graph ) && ((connecting_val >= 2) || 
                                  (connecting_val >= 1  && 
                                   isaltsplice(motif)))
          psi = Nullable( 0.99 ) #&& likelihood_ci( psi, inc_cnt, z=1.64 )
      else
         # NA is ignored.
      end
   else # there is skipping
      if isnull( inc_graph )
         if spanning_val >= 1
             psi = Nullable( 0.0 ) #&& likelihood_ci( psi, exc_graph.count, z=1.64 )
         else
            # NA
         end
      else
         get(exc_graph).psi = zeros( length( get(exc_graph).count ) )
         get(inc_graph).psi = zeros( length( get(inc_graph).count ) )
         calculate_psi!( inc_graph.value, exc_graph.value, [get(inc_graph).count; get(exc_graph).count] )
         #println( "$inc_path\n$exc_graph\n$ambig_cnt\n\n" )
         it = rec_spliced_em!( inc_graph.value, exc_graph.value, ambig_cnt.value, sig=4 )
         #println( "$inc_path\n$exc_graph\n$ambig_cnt\n\n" )
         psi = Nullable( sum( get(inc_graph).psi ) )
      end
   end
   psi,inc_graph,exc_graph,ambig_cnt
end

function extend_edges!{K,V}( edges::IntervalMap{K,V}, pgraph::PsiGraph, igraph::PsiGraph,
                             ambig_edge::Nullable{Vector{IntervalValue}}, node::NodeInt )

   minv = min( pgraph.min, igraph.min )
   maxv = max( pgraph.max, igraph.max  )
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
            if edg.last in igraph
               push!( igraph, edg, value_bool=false )
               shouldpush = true
            end
            if shouldpush
               if edg.first < minv
                  minv = edg.first
               end
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
            if edg.first in pgraph
               push!( pgraph, edg, value_bool=false )
               shouldpush = true
            end
            if edg.first in igraph
               push!( igraph, edg, value_bool=false )
               shouldpush = true
            end
            if shouldpush
               if edg.last > maxv
                  maxv = edg.last
               end
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

function process_events( outfile, lib::GraphLib, graphq::GraphLibQuant; isnodeok=true )
   io = open( outfile, "w" )
   stream = ZlibDeflateOutputStream(io)
   for g in 1:length(lib.graphs)
      name = lib.names[g]
      chr  = lib.info[g].name
      strand = lib.info[g].strand ? '+' : '-'
      #println(STDERR, "$g, $name, $chr, $strand" )
      _process_events( stream, lib.graphs[g], graphq.quant[g], (name,chr,strand), isnodeok=isnodeok )
   end
   close(stream)
   close(io)
end

function _process_events( io::BufOut, sg::SpliceGraph, sgquant::SpliceGraphQuant, info::GeneMeta; isnodeok=false )
   # Heres the plan:
   # step through sets of edges, look for edge motifs, some are obligate calculations
   # others we only calculate psi if there is an alternative edge
   # Then calculate bias of the event 
   if length( sg.edgetype ) <= 2
      return 0
   end
   i = 1
   while i < length(sg.edgetype)
      motif = convert(EdgeMotif, sg.edgetype[i], sg.edgetype[i+1] )
      motif == NONE_MOTIF && (i += 1; continue)
      if isobligate( motif ) # is utr event
         psi,utr,ambig = _process_tandem_utr( sg, sgquant, convert(NodeInt, i), motif ) 
         if !isnull( psi )
            if any( map( isnan, psi.value ) )
               println(STDERR, get(utr))
            end
            ambig_cnt = isnull( ambig ) ? 0.0 : sum( ambig.value )
            i = output_utr( io, round(get(psi),4), utr, ambig_cnt, motif, sg, i , info )
            #i += (motif == TXST_MOTIF) ? length(psi.value)-2 : length(psi.value)-2 
         end
      else  # is a spliced node
         bias = calculate_bias!( sgquant )
         psi,inc,exc,ambig = _process_spliced( sg, sgquant, convert(NodeInt, i), motif, bias, isnodeok )
         if !isnull( psi )
            total_cnt = sum(inc) + sum(exc) + sum(ambig)
            conf_int  = binomial_likelihood_ci( get(psi), total_cnt, sig=3 )
            output_psi( io, signif(get(psi),4), inc, exc, total_cnt, conf_int, motif, sg, i, info, bias  ) # TODO bias
         end
      end
      i += 1
   end
   # process back-splicing
   output_circular( io, sg, sgquant, info )
end

function Base.unsafe_copy!{T <: Number}( dest::Vector{T}, src::Vector{T}; indx_shift=0 )
   for i in 1:length(src)
      dest[i+indx_shift] = src[i]
   end
end

function divsignif!{ N <: Number, D <: Number, I <: Integer }( arr::Vector{N}, divisor::D, sig::I )
   if sig > 0
      for i in 1:length( arr )
         @fastmath arr[i] = signif( arr[i] / divisor, sig )
      end
   else
      for i in 1:length( arr )
         @fastmath arr[i] /= divisor
      end
   end
end

@inline function binomial_likelihood_ci( p, n, z=1.64; sig=0 )
   const fisher_info = (p * (1-p)) / n
   if fisher_info < 0
      return(0.0,1.0)
   end
   const ci = z * sqrt( fisher_info )
   if sig > 0
      const lo = signif( max( 0.0, p - ci ), sig )
      const hi = signif( min( 1.0, p + ci ), sig )
   else
      const lo = max( 0.0, p - ci )
      const hi = min( 1.0, p + ci )
   end
   lo,hi
end

# spliced psi
function calculate_psi!( igraph::PsiGraph, egraph::PsiGraph, counts::Vector{Float64}; sig=0 )
   for i in 1:length(igraph.psi)
      igraph.psi[i] = counts[i] / igraph.length[i]
   end
   for i in 1:length(egraph.psi)
      idx = i + length(igraph.psi)
      egraph.psi[i] = counts[idx] / egraph.length[i]
   end
   cnt_sum = sum( egraph.psi ) + sum( igraph.psi )
   divsignif!( igraph.psi, cnt_sum, sig )
   divsignif!( egraph.psi, cnt_sum, sig )
end

function calculate_psi!( pgraph::PsiGraph, counts::Vector{Float64}; sig=0 )
   for i in 1:length(pgraph.psi)
      pgraph.psi[i] = counts[i] / pgraph.length[i]
   end
   const cnt_sum = sum( pgraph.psi )
   divsignif!( pgraph.psi, cnt_sum, sig )
end

function rec_spliced_em!( igraph::PsiGraph, egraph::PsiGraph, 
                          ambig::Vector{AmbigCounts};
                          inc_temp::Vector{Float64}=zeros(length(igraph.count)),
                          exc_temp::Vector{Float64}=zeros(length(egraph.count)),
                          count_temp::Vector{Float64}=ones(length(igraph.count)+length(egraph.count)),
                          it=1, max=1000, sig=0 )

   unsafe_copy!( count_temp, igraph.count, indx_shift=0 )
   unsafe_copy!( count_temp, egraph.count, indx_shift=length(igraph.count) )
   unsafe_copy!( inc_temp, igraph.psi )
   unsafe_copy!( exc_temp, egraph.psi )
   
   for ac in ambig
      idx = 1
      if it > 1 # maximization
         ac.prop_sum = 0.0
         for p in ac.paths
            prev_psi = p <= length(igraph.psi) ? igraph.psi[p] : egraph.psi[p-length(igraph.psi)]
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

   calculate_psi!( igraph, egraph, count_temp, sig=sig ) # expectation

   if (inc_temp != igraph.psi || egraph.psi != exc_temp) && it < max
      it = rec_spliced_em!( igraph, egraph, ambig, 
                            inc_temp=inc_temp, exc_temp=exc_temp, count_temp=count_temp,
                            it=it+1, max=max, sig=sig )
   end
   it
end

function rec_tandem_em!( pgraph::PsiGraph, ambig::Vector{AmbigCounts};
                          utr_temp::Vector{Float64}=zeros(length(pgraph.count)),
                          count_temp::Vector{Float64}=ones(length(pgraph.count)),
                          it=1, max=1000, sig=0 )

   unsafe_copy!( count_temp, pgraph.count, indx_shift=0 )
   unsafe_copy!( utr_temp, pgraph.psi )

   for ac in ambig
      idx = 1
      if it > 1 # maximization
         ac.prop_sum = 0.0
         for p in ac.paths
            prev_psi = pgraph.psi[p] 
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

   calculate_psi!( pgraph, count_temp, sig=sig ) # expectation

   if utr_temp != pgraph.psi && it < max
      it = rec_tandem_em!( pgraph, ambig,
                            utr_temp=utr_temp, count_temp=count_temp,
                            it=it+1, max=max, sig=sig )
   end
   it
end


