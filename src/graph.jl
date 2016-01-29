
using IntervalTrees

typealias Exonmax UInt16
typealias Edgetype UInt8 # Encode Edgetype as:
# 'L' = 0x05
# 'R' = 0x06
# 'S' = 0x07, therfore... 'SL' = 0x07-0x05 = 0x02, 'SR' = 0x01

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
  seq::DNASequence
end

# empty constructor
SpliceGraph() = SpliceGraph( Vector{Coordint}(), Vector{Coordint}(),
                             Vector{Exonmax}(),  Vector{Exonmax}(),
                             Vector{Exonmax}() )

# Main constructor
# Build splice graph here.
function SpliceGraph( gene::Refgene, chrom::DNASequence )
   # splice graph variables
   nodeoffset = Vector{Coordint}()
   nodesize   = Vector{Coordint}()
   edgeoffset = Vector{Coordint}()
   edgetype   = Vector{Edgetype}()
   seq        = DNASequence()
   

   for i in 1:length(don)
      # iterate through donors, and acceptors
      # left to right. '-' = rc unshift?
   end
   return SpliceGraph( nodeoffset, edgeoffset, 
                       left, right, txend, seq )
end
