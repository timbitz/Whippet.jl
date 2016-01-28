
using IntervalTrees

typealias ExonMax UInt16

type SpliceGraph
   nodecnt::Vector{Float64}
   edgecnt::IntervalMap{ExonMax,Float64}
   circcnt::Dict{Tuple{ExonMax,ExonMax},Float64}
end

SpliceGraph() = SpliceGraph(Vector{Float64}(), 
                            IntervalMap{ExonMax,Float64}(),
                            Dict{Tuple{ExonMax,ExonMax},Float64}())


