# Extended from events.jl for full isoforms

function full_graph( sgquant::SpliceGraphQuant )
   graph = Nullable{PsiGraph}()
   for edg in sgquant.edge
      if isnull( graph )
         graph = Nullable(PsiGraph( Vector{Float64}(), Vector{Float64}(),
                                    Vector{Float64}(), Vector{IntSet}(),
                                    edg.first, edg.last ))
      end
      push!( graph.value, edg, value_bool=false )
   end
   graph
end

function process_isoforms( sg::SpliceGraph, sgquant::SpliceGraphQuant )
   graph = full_graph( sgquant )
   if !isnull( graph )
      graph = Nullable( reduce_graph( graph.value ) )
   end
    
end

# TODO: filter isoforms with out txStart and End
# TODO: find largest ORF

# SpliceGraph now records annotated edges in an ExonTree where the intervals
# are encoded as (1,2) for an annotated edge between nodes 1 and 2
# returns: nothing
function add_path_edges!( edges::ExonTree, path::IntSet )
   s = start(path)
   lastv,s = next( path, s )
   while !done( path, s )
      nextv,s = next( path, s )
      const interv = Interval{ExonInt}( lastv, nextv )
      push!( edges, interv )
      lastv = nextv
   end
end

# Build minimal set of paths to explain graph
# Explanation: If we have 3 edges for example 1-2, 2-3, and 1-3
# We have to create 2 paths, 1-2-3, and 1-3, since 1-3 can't
# possibly contain the second node.
# This function will do this with any number of edges or complex
# paths through a graph and is the essence of how Whippet quantification
# works.. these graph paths are then assigned unambiguous counts, while
# ambiguous counts are then assigned to each graph path using the EM algorithm
function reduce_graph( edges::ExonTree )
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
