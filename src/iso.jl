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
# TODO: 
