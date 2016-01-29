
using IntervalTrees

typealias Exonmax UInt16

bitstype 8 EdgeType
const EDGETYPE_TO_UINT8 = fill( 0x07, (8,8) )
EDGETYPE_TO_UINT8[8,6] = 0x00 # 'SL' = 0x00 # Tx Start
EDGETYPE_TO_UINT8[8,7] = 0x01 # 'SR' = 0x01 # Tandem Last Exon
EDGETYPE_TO_UINT8[6,8] = 0x02 # 'LS' = 0x02 # Tandem First Exon
EDGETYPE_TO_UINT8[7,8] = 0x03 # 'RS' = 0x03 # Tx End
EDGETYPE_TO_UINT8[6,7] = 0x04 # 'LR' = 0x04 # Intron 5'->3' SS
EDGETYPE_TO_UINT8[6,6] = 0x05 # 'LL' = 0x05 # Alt- 5'SS
EDGETYPE_TO_UINT8[7,7] = 0x06 # 'RR' = 0x06 # Alt- 3'SS

function Base.convert( ::Type{EdgeType}, one::UInt8, two::UInt8 )
   @assert( 5 <= one <= 7 && 5 <= one <= 7 ) 
   EDGETYPE_TO_UINT8[one+1,two+1]
end

Base.convert( ::Type{EdgeType}, one::DNANucleotide, two::DNANucleotide ) = Base.convert( EdgeType, convert(UInt8, one), convert(UInt8, two) )
Base.convert( ::Type{EdgeType}, edge::UInt8 ) = box(EdgeType, unbox(UInt8, edge ))
Base.convert( ::Type{UInt8}, edge::EdgeType ) = box(UInt8, unbox(EdgeType, edge ))

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
  nodelen::Vector{Coordint}
  edgetype::Vector{EdgeType}
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
   nodelen    = Vector{Coordint}()
   edgetype   = Vector{EdgeType}()
   seq        = DNASequence()
   
   a = length(gene.acc)
   d = length(gene.don)

   for i in 1:length(don)
      # iterate through donors, and acceptors
      # left to right. '-' = rc unshift?
   end
   return SpliceGraph( nodeoffset, nodelen, edgetype, seq )
end
