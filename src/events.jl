
const COMPLEX_CHAR = 'K'

primitive type EdgeMotif 8 end

Base.convert( ::Type{EdgeMotif}, motif::UInt8 ) = reinterpret(EdgeMotif, motif )
# Base.convert( ::Type{UInt8}, motif::EdgeMotif ) = reinterpret(UInt8, motif )

UInt8(motif::EdgeMotif) = reinterpret(UInt8, motif)

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
      MOTIF_STRING[ UInt8(NONE_MOTIF) + 1 ] = "NA"

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

Base.convert(::Type{S}, edg::EdgeMotif ) where {S <: AbstractString} = MOTIF_STRING[ UInt8(edg) + 1 ]

#isobligate(  motif::EdgeMotif ) = motif != NONE_MOTIF && !( UInt8(motif) & 0b100 == 0b100 )
isobligate( motif::EdgeMotif ) = (0 <= UInt8(motif) <= 1)
isaltsplice( motif::EdgeMotif ) = (UInt8(motif) & 0b110) == 0b110

isspanning( edge::I, node::T ) where {I <: AbstractInterval, T <: Integer} = edge.first < node < edge.last ? true : false
isconnecting( edge::I, node::T ) where {I <: AbstractInterval, T <: Integer} = edge.first == node || edge.last == node ? true : false

# this is meant for short arrays when it is faster
# than using the overhead of a set
# DEPRECATED, BitSet is very efficient
function unique_push!( arr::Vector{T}, el::T ) where T
   if !( el in arr )
      push!( arr, el )
   end
end


# This holds one or many sets of connected
# nodes + the count of the reads and the eff_len
# We also store the min and max node of all sets
# Note: It did occur to me that a cleaner solution
# is PsiGraph holding paths::Vector{PsiPath} (deprecated)
# but the current implementation was choosen for
# faster memory access w.r.t. stride.
mutable struct PsiGraph
   psi::Vector{Float64}
   count::Vector{Float64}
   length::Vector{Float64}
   nodes::Vector{BitSet}
   min::NodeInt
   max::NodeInt
end

Base.sum( pgraph::PsiGraph ) = sum( pgraph.count )
Base.sum( ngraph::Nullable{PsiGraph} ) = isnull( ngraph ) ?
                                         0.0 : convert(Float64, sum( ngraph.value ))

function hasintersect( a::BitSet, b::BitSet )
   seta,setb = length(a) > length(b) ? (b,a) : (a,b)
   for elem in seta
      if elem in setb
         return true
      end
   end
   false
end

function hasintersect_terminal( a::BitSet, b::BitSet )
   if first(a) == last(b) || first(b) == last(a)
      return true
   else
      return false
   end
end

function hasintersect_terminal( a, b )
   if first(a) == last(b) || first(b) == last(a)
      return true
   else
      return false
   end
end

function Base.in( i::T, pgraph::PsiGraph ) where T <: Integer
   for nodeset in pgraph.nodes
      (i in nodeset) && return true
   end
   return false
end

complexity( one::PsiGraph, two::PsiGraph ) = complexity(length(one.nodes) + length(two.nodes))
complexity( one::PsiGraph ) = complexity(length(one.nodes))
# This function calculates complexity of a splicing event
# as the ceil log2 number of paths through the graph
complexity( num_paths::Int ) = @fastmath Int(ceil(log2( num_paths )))

# Calculate Shannon's Entropy, -Sum( Pi * log2( Pi ) )
function shannon_index( probs::Vector{Float64} )
   index = 0.0
   for i in 1:length(probs)
      prob = (probs[i] == 0.0 ? eps(Float64) : probs[i])
      index += prob * log2(prob)
   end
   index * -1
end

shannon_index( one::PsiGraph ) = shannon_index( one.psi )

# The Psi vectors in these two PsiGraphs need to sum to one!
# this does not assert this!
function shannon_index( one::PsiGraph, two::PsiGraph )
   index = 0.0
   for i in 1:length(one.psi)
      prob = (one.psi[i] == 0.0 ? eps(Float64) : one.psi[i])
      index += prob * log2(prob)
   end
   for j in 1:length(two.psi)
      prob = (two.psi[j] == 0.0 ? eps(Float64) : two.psi[j])
      index += prob * log2(prob)
   end
   index * -1
end


# Build minimal set of paths to explain graph
# Explanation: If we have 3 edges for example 1-2, 2-3, and 1-3
# We have to create 2 paths, 1-2-3, and 1-3, since 1-3 can't
# possibly contain the second node.
# This function will do this with any number of edges or complex
# paths through a graph and is the essence of how Whippet quantification
# works.. these graph paths are then assigned unambiguous counts, while
# ambiguous counts are then assigned to each graph path using the EM algorithm

function reduce_graph( edges::PsiGraph, graph::PsiGraph=deepcopy(edges), maxit::Int=100, maxiso::Int=1024 )
   newgraph = PsiGraph( Vector{BitSet}(),  Vector{Float64}(),
                        Vector{Float64}(), Vector{Float64}(),
                        graph.min, graph.max )
   it = 1
   while (newgraph.nodes != graph.nodes || it == 1) && it <= maxit
      if it > 1
         graph,newgraph = newgraph,graph
         empty!( newgraph.nodes  )
         empty!( newgraph.length )
         empty!( newgraph.count  )
         if length(graph.nodes) > maxiso
            ord = shuffle(1:length(graph.nodes))
            graph.nodes  = graph.nodes[ord[1:maxiso]]
            graph.length = graph.length[ord[1:maxiso]]
            graph.count  = graph.count[ord[1:maxiso]]
         end
      end
      @inbounds for i in 1:length(graph.nodes)
         push_i = false
         for j in 1:length(edges.nodes)
            if hasintersect_terminal( graph.nodes[i], edges.nodes[j] )
               push_i = true
               jointset = union( graph.nodes[i], edges.nodes[j] )
               (jointset in newgraph.nodes) && continue
               push!( newgraph.nodes,  jointset )
               push!( newgraph.length, graph.length[i] + edges.length[j] )
               push!( newgraph.count,  graph.count[i] + edges.count[j] )
               newgraph.min = min( min(first(graph.nodes[i]), first(edges.nodes[j])), newgraph.min )
               newgraph.max = max( max(last(graph.nodes[i]),  last(edges.nodes[j])),  newgraph.max )
            end
         end
         if !push_i && !(graph.nodes[i] in newgraph.nodes)
            push!( newgraph.nodes, graph.nodes[i] )
            push!( newgraph.length, graph.length[i] )
            push!( newgraph.count, graph.count[i] )
            newgraph.min = min( first(graph.nodes[i]), newgraph.min )
            newgraph.max = max( last(graph.nodes[i]),  newgraph.max )
         end
      end
      it += 1
   end
   newgraph
end

@inline function reduce_graph_simple( edges::PsiGraph )
   if length(edges.nodes) == 2
      if hasintersect_terminal( edges.nodes[1], edges.nodes[2] )
         union!( edges.nodes[1], edges.nodes[2] )
         edges.length[1] += edges.length[2]
         edges.count[1]  += edges.count[2]
         pop!( edges.nodes )
         pop!( edges.length )
         pop!( edges.count )
      end
      return edges
   else
      return reduce_graph( edges )
   end
end


function Base.push!( pgraph::PsiGraph, edg::I;
                     value_bool=true, length=1.0 ) where I <: AbstractInterval
   push!( pgraph.count, value_bool ? get(edg.value) : 0.0 )
   push!( pgraph.length, value_bool ? length : 0.0 )
   push!( pgraph.nodes, BitSet([edg.first, edg.last]) )
   #reduce_graph_terminal!( pgraph )
   if edg.first < pgraph.min
      pgraph.min = edg.first
   end
   if edg.last > pgraph.max
      pgraph.max = edg.last
   end
end


mutable struct AmbigCounts
   paths::Vector{NodeInt}
   prop::Vector{Float64}
   prop_sum::Float64
   multiplier::Float64
end

Base.:(==)( a::AmbigCounts, b::AmbigCounts ) = a.paths == b.paths
Base.:(==)( a::AmbigCounts, b::BitSet ) = ( a.paths == b )
Base.:(==)( a::BitSet, b::AmbigCounts ) = ( a == b.paths )
Base.:(==)( a::Vector{I}, b::BitSet ) where {I <: Integer} = ( b == a )
function Base.:(==)( iset::BitSet, ivec::Vector{I} ) where I <: Integer
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

function extend_edges!( edges::IntervalMap{K,V}, pgraph::PsiGraph, igraph::PsiGraph,
                        ambig_edge::Nullable{Vector{IntervalValue}}, node::NodeInt ) where {K,V}

   minv = min( pgraph.min, igraph.min )
   maxv = max( pgraph.max, igraph.max )
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

# add_node_counts! for single pgraph (used in tandem_utrs)
function add_node_counts!( ambig::Vector{AmbigCounts}, pgraph::PsiGraph,
                           sgquant::SpliceGraphQuant, bias::Float64 )

   iset = BitSet()
   for n in pgraph.min:pgraph.max # lets go through all possible nodes
      for i in 1:length(pgraph.nodes)
         if n in pgraph.nodes[i]
            push!( iset, i )
            pgraph.length[i] += sgquant.leng[n]
         end
      end

      #=                                            =#
      length( iset ) == 0 && continue # unspliced node
      if length( iset ) == 1 # non-ambiguous node
         if first( iset ) <= length(pgraph.nodes) # belongs to inclusion
            idx = first( iset )
            @fastmath pgraph.count[idx] += get(sgquant.node[n]) #* bias / sgquant.leng[n]
         end
      else
         #check if there is already an entry for this set of paths
         # if so, just increment multiplier, if not, make new one
         exists = false
         for am in ambig
            if iset == am
               @fastmath am.multiplier += get(sgquant.node[n]) #* bias / sgquant.leng[n]
               exists = true
               break
            end
         end
         if !exists && get(sgquant.node[n]) > 0
            push!( ambig, AmbigCounts( collect(iset),
                                       ones( length(iset) ) / length(iset),
                                       1.0, get(sgquant.node[n]) ) ) #* bias / sgquant.leng[n] ) )
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
   iset = BitSet()
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
            @fastmath igraph.count[idx] += get(sgquant.node[n]) * bias / sgquant.leng[n]
         else # belongs to exclusion
            idx = first( iset ) - length(igraph.nodes)
            @fastmath egraph.count[idx] += get(sgquant.node[n]) * bias / sgquant.leng[n]
         end
      else
         #check if there is already an entry for this set of paths
         # if so, just increment multiplier, if not, make new one
         exists = false
         for am in ambig
            if iset == am
               @fastmath am.multiplier += get(sgquant.node[n]) * bias / sgquant.leng[n]
               exists = true
               break
            end
         end
         if !exists && get(sgquant.node[n]) > 0
            push!( ambig, AmbigCounts( collect(iset),
                                       ones( length(iset) ) / length(iset),
                                       1.0, get(sgquant.node[n]) * bias / sgquant.leng[n] ) )
         end
      end
      empty!( iset ) # clean up
   end
end

# add_edge_counts! for one psigraph, (used in tandem utr)
function add_edge_counts!( ambig::Vector{AmbigCounts}, pgraph::PsiGraph,
                           sgquant::SpliceGraphQuant )
   iset = BitSet()
   for edg in sgquant.edge
      for i in 1:length(pgraph.nodes)
         if edg in pgraph.nodes[i]
            push!(iset, i)
            #pgraph.length[i] += 1
         end
      end
      #=                                             =#
      length( iset ) == 0 && continue
      if length( iset ) == 1
         idx = first( iset )
         pgraph.count[idx] += get(edg.value)
      else
         exists = false
         for am in ambig
            if iset == am
               am.multiplier += get(edg.value)
               exists = true
               break
            end
         end
         if !exists && get(edg.value) > 0
            push!( ambig, AmbigCounts( collect(iset),
                                       ones( length(iset) ) / length(iset),
                                       1.0, get(edg.value) ) )
         end
      end
      empty!( iset )
   end
end

# add_edge_counts! for two graphs, inc_graph and exc_graph! (used in spliced events)
function add_edge_counts!( ambig::Vector{AmbigCounts}, igraph::PsiGraph,
                           egraph::PsiGraph, edges::Vector{IntervalValue} )
   iset = BitSet()
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
            igraph.count[idx] += get(edg.value)
         else
            idx = first( iset ) - length(igraph.nodes)
            egraph.count[idx] += get(edg.value)
         end
      else
         exists = false
         for am in ambig
            if iset == am
               am.multiplier += get(edg.value)
               exists = true
               break
            end
         end
         if !exists && get(edg.value) > 0
            push!( ambig, AmbigCounts( collect(iset),
                                       ones( length(iset) ) / length(iset),
                                       1.0, get(edg.value) ) )
         end
      end
      empty!( iset )
   end
end

function build_utr_graph( nodes::BitSet, motif::EdgeMotif, sgquant::SpliceGraphQuant )

   utr_graph = Nullable(PsiGraph( Vector{Float64}(), Vector{Float64}(),
                                  Vector{Float64}(), Vector{BitSet}(),
                                  first(nodes), last(nodes) ))
   curset = BitSet()

   while length(nodes) > 0
      if motif == TXST_MOTIF
         curval = pop!(nodes)
      else
         curval = popfirst!(nodes)
      end
      push!( curset, curval )
      push!( get(utr_graph).nodes, copy(curset) )
      push!( get(utr_graph).psi, 0.0 )
      push!( get(utr_graph).length, 0.0 )
      push!( get(utr_graph).count, 0.0 )
   end

   utr_graph
end

function _process_tandem_utr( sg::SpliceGraph, sgquant::SpliceGraphQuant{C,R},
                              node::NodeInt, motif::EdgeMotif, readlen::Int ) where {C <: SGAlignContainer, R <: ReadCounter}

   utr_graph  = Nullable{PsiGraph}()
   ambig_cnt  = Nullable{Vector{AmbigCounts}}()

   used_node  = BitSet()
   total_cnt  = 0.0

   if motif == TXEN_MOTIF
      # add node directly upstream (CE next to tandemUTR)
      push!( used_node, node-1 )
      total_cnt += get(sgquant.node[node-1])
   end

   i = node
   curmotif = motif
   while curmotif == motif
      push!( used_node, i )
      total_cnt += get(sgquant.node[i])
      if i > 1
         interv = IntervalTrees.Interval{ExonInt}( i-1, i )
         total_cnt += get(get( sgquant.edge, interv, IntervalValue(0,0,default(R)) ).value)
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
      total_cnt += get(sgquant.node[i])
      interv = IntervalTrees.Interval{ExonInt}( i-1, i )
      total_cnt += get(get( sgquant.edge, interv, IntervalValue(0,0,default(R)) ).value)
   end

   psi = Nullable{Float64}()
   len = length(used_node)
   if total_cnt > 2
      ambig_cnt = Nullable( Vector{AmbigCounts}() )
      utr_graph = build_utr_graph( used_node, motif, sgquant )
      add_node_counts!( ambig_cnt.value, utr_graph.value, sgquant, 1.0 )
      add_edge_counts!( ambig_cnt.value, utr_graph.value, sgquant )
      it = rec_tandem_em!( utr_graph.value, ambig_cnt.value, sig=4, readlen=readlen )
      psi = Nullable( get(utr_graph).psi )
   end

   psi,utr_graph,ambig_cnt,len
end

function _process_spliced( sg::SpliceGraph, sgquant::SpliceGraphQuant,
                           node::NodeInt, motif::EdgeMotif, bias::Float64, isnodeok::Bool,
                           minedgeweight::Float64=0.02 )

   inc_graph  = Nullable{PsiGraph}()
   exc_graph  = Nullable{PsiGraph}()
   ambig_edge = Nullable{Vector{IntervalValue}}()
   ambig_cnt  = Nullable{Vector{AmbigCounts}}()

   connecting_val = 0.0
   spanning_val   = 0.0

   maxvalue = length(sgquant.edge) > 0 ? get(maximum( sgquant.edge ).value) : 0.0

   for edg in intersect( sgquant.edge, (node, node) )

      (get(edg.value) / maxvalue < minedgeweight) && continue

      if   isconnecting( edg, node )
         if isnull( inc_graph )
            inc_graph = Nullable(PsiGraph( Vector{Float64}(), Vector{Float64}(),
                                           Vector{Float64}(), Vector{BitSet}(),
                                           edg.first, edg.last ))
         end
         if isnull( ambig_edge )
            ambig_edge = Nullable( Vector{IntervalValue}() )
         end
         push!( inc_graph.value, edg, value_bool=false )
         push!( ambig_edge.value, edg )
         connecting_val += get(edg.value)
      elseif isspanning( edg, node )
         if isnull( exc_graph ) #don't allocate unless there is alt splicing
            exc_graph = Nullable(PsiGraph( Vector{Float64}(), Vector{Float64}(),
                                           Vector{Float64}(), Vector{BitSet}(),
                                           edg.first, edg.last ))
         end
         if isnull( ambig_edge )
            ambig_edge = Nullable( Vector{IntervalValue}() )
         end
         push!( exc_graph.value, edg, value_bool=false )
         push!( ambig_edge.value, edg )
         spanning_val += get(edg.value)
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
   else
      !isnull( inc_graph ) && (inc_graph = Nullable( reduce_graph_simple( inc_graph.value ) ))
      !isnull( exc_graph ) && (exc_graph = Nullable( reduce_graph_simple( exc_graph.value ) ))
   end # end expanding module

  # lets finish up now.
   psi = Nullable{Float64}()
   if isnull( exc_graph ) # no spanning edge
      # check if we have both inclusion edges represented, or one if alt 5'/3'
      if !isnull( inc_graph ) && ((connecting_val >= 2) ||
                                  (connecting_val >= 1  &&
                                   isaltsplice(motif)))
          psi = Nullable( 1.0 ) #&& likelihood_ci( psi, inc_cnt, z=1.64 )
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
         it = spliced_em!( inc_graph.value, exc_graph.value, ambig_cnt.value, sig=4 )
         psi = Nullable( sum( get(inc_graph).psi ) )
      end
   end
   total_reads = connecting_val + spanning_val
   psi,inc_graph,exc_graph,ambig_cnt,total_reads
end

function process_events( outfile, lib::GraphLib, graphq::GraphLibQuant;
                         isnodeok=false, iscircok=false, readlen::Int=50 )
   io = open( outfile, "w" )
   stream = ZlibDeflateOutputStream( io )
   output_psi_header( stream )
   for g in 1:length(lib.graphs)
      name = lib.names[g]
      chr  = lib.info[g].name
      strand = lib.info[g].strand ? '+' : '-'
      _process_events( stream, lib.graphs[g], graphq.quant[g], (name,chr,strand),
                       isnodeok=isnodeok, iscircok=iscircok, readlen=readlen )
   end
   close(stream)
   close(io)
end

function _process_events( io::BufOut, sg::SpliceGraph, sgquant::SpliceGraphQuant, info::GeneMeta;
                          isnodeok=false, iscircok=false, readlen::Int=50.0 )
   # Heres the plan:
   # step through sets of edges, look for edge motifs, some are obligate calculations
   # others we only calculate psi if there is an alternative edge
   if length( sg.edgetype ) <= 2
      return 0
   end
   i = 1
   while i < length(sg.edgetype)
      motif = convert(EdgeMotif, sg.edgetype[i], sg.edgetype[i+1] )
      next_istandem = i+1 < length(sg.edgetype) ? isobligate( convert(EdgeMotif, sg.edgetype[i+1], sg.edgetype[i+2]) ) : false
      if motif == NONE_MOTIF
         if !next_istandem
            output_empty( io, motif, sg, i, info )
         end
      elseif isobligate( motif ) # is utr event
         psi,utr,ambig,len = _process_tandem_utr( sg, sgquant, convert(NodeInt, i), motif, readlen )
         if !isnull( psi ) && !any( map( isnan, psi.value ) )
            total_cnt = sum(utr) + sum(ambig)
            i = output_utr( io, round.(get(psi), digits=4), utr, total_cnt, motif, sg, i , info )
         else
            # psi/utr/total_cnt ignored here.
            i = output_utr( io, zeros(len), utr, 0.0, motif, sg, i, info, empty=true )
         end
      else  # is a spliced node
         bias = calculate_bias!( sgquant )
         psi,inc,exc,ambig,total_cnt = _process_spliced( sg, sgquant, convert(NodeInt, i), motif, bias, isnodeok )
         #total_cnt = sum(inc) + sum(exc) + sum(ambig)
         if !isnull( psi ) && 0 <= psi.value <= 1 && total_cnt > 0
            conf_int  = beta_posterior_ci( psi.value, total_cnt, sig=3 )
            output_psi( io, round(psi.value, sigdigits=3), inc, exc, total_cnt, conf_int, motif, sg, sgquant.edge, i, info, bias  ) # TODO bias
         else
            output_empty( io, motif, sg, i, info )
         end
      end
      i += 1
   end
   # process back-splicing
   iscircok && output_circular( io, sg, sgquant, info )
end

function divsignif!( arr::Vector{N}, divisor::D, sig::I ) where { N <: Number, D <: Number, I <: Integer }
   if sig > 0
      for i in 1:length( arr )
         @fastmath arr[i] = round( arr[i] / divisor, sigdigits=sig )
      end
   else
      for i in 1:length( arr )
         @fastmath arr[i] /= divisor
      end
   end
end


# Beta(
@inline function beta_posterior_ci( p, n; ci=0.9, sig=0 )
   lo_q = (1 - ci)/2
   hi_q = 1 - lo_q
   beta = Beta( p*n + 1, (1-p)*n + 1 )
   if sig > 0
      lo = round( quantile(beta, lo_q), sigdigits=sig )
      hi = round( quantile(beta, hi_q), sigdigits=sig )
   else
      lo = quantile(beta, lo_q)
      hi = quantile(beta, hi_q)
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

function calculate_psi!( pgraph::PsiGraph, counts::Vector{Float64}; sig=0, readlen::Int=50)
   for i in 1:length(pgraph.psi)
      pgraph.psi[i] = counts[i] / max( pgraph.length[i] - readlen, 1 )
   end
   cnt_sum = sum( pgraph.psi )
   divsignif!( pgraph.psi, cnt_sum, sig )
end


function spliced_em!( igraph::PsiGraph, egraph::PsiGraph, ambig::Vector{AmbigCounts};
                      it=1, maxit=1500, sig=0 )

   inc_temp   = zeros(length(igraph.count))
   exc_temp   = zeros(length(egraph.count))
   count_temp = ones(length(igraph.count)+length(egraph.count))

   while (inc_temp != igraph.psi || egraph.psi != exc_temp) && it < maxit

      unsafe_copy!( count_temp, igraph.count, indx_shift=0 )
      unsafe_copy!( count_temp, egraph.count, indx_shift=length(igraph.count) )
      unsafe_copy!( inc_temp, igraph.psi )
      unsafe_copy!( exc_temp, egraph.psi )

      for ac in ambig
         idx = 1
         if it > 1
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

      calculate_psi!( igraph, egraph, count_temp, sig=sig )

      it += 1
   end
   it
end

function rec_tandem_em!( pgraph::PsiGraph, ambig::Vector{AmbigCounts};
                          utr_temp::Vector{Float64}=zeros(length(pgraph.count)),
                          count_temp::Vector{Float64}=ones(length(pgraph.count)),
                          it=1, maxit=1500, sig=0, readlen::Int=50 )

   unsafe_copy!( count_temp, pgraph.count, indx_shift=0 )
   unsafe_copy!( utr_temp, pgraph.psi )

   for ac in ambig
      idx = 1
      if it > 1
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

   calculate_psi!( pgraph, count_temp, sig=sig, readlen=readlen )

   if utr_temp != pgraph.psi && it < maxit
      it = rec_tandem_em!( pgraph, ambig,
                           utr_temp=utr_temp, count_temp=count_temp,
                           it=it+1, maxit=maxit, sig=sig, readlen=readlen )
   end
   it
end
