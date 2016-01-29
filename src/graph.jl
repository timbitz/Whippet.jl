
using IntervalTrees

typealias Exonmax UInt16

# This is where we count reads for nodes/edges/circular-edges
type SpliceGraphQuant 
   nodecnt::Vector{Float64}
   edgecnt::IntervalMap{Exonmax,Float64}
   circcnt::Dict{Tuple{Exonmax,Exonmax},Float64}
end

# Default constructer
SpliceGraphQuant() = SpliceGraphQuant( Vector{Float64}(), Vector{UInt32}(),
                                       IntervalMap{Exonmax,Float64}(),
                                       Dict{Tuple{Exonmax,Exonmax},Float64}() )


# This holds a representation of the splice graph
# which is a directed multigraph
immutable SpliceGraph
  nodeoffset::Vector{Coordint}
  edgeoffset::Vector{Coordint}
  left::Vector{Exonmax} #donor splice sites indexes
  right::Vector{Exonmax} #acceptor splice site indexes
  txend::Vector{Exonmax}  #tx starts/ends
end

SpliceGraph() = SpliceGraph( Vector{Coordint}(), Vector{Coordint}(),
                             Vector{Exonmax}(), Vector{Exonmax}(),
                             Vector{Exonmax}() )

function SpliceGraph( gene::Refgene )
   ss = SpliceGraph()
   for d in gene.don

   end  
end
