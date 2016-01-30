
# requires:
# include("bio_nuc_patch.jl")
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
  nodeoffset::Vector{Coordint} # SG offset
  nodecoord::Vector{Coordint}  # Genome offset
  nodelen::Vector{Coordint}
  edgetype::Vector{EdgeType}
  seq::DNASequence
end
# All positive strand oriented sequences---> 
# Node array: txStart| 1 |   2   | 3 |    4    |5| 6 |txEnd
# Edge array:        1   2       3   4         5 6   7

# empty constructor
SpliceGraph() = SpliceGraph( Vector{Coordint}(), Vector{Coordint}(),
                             Vector{Exonmax}(),  Vector{Exonmax}(),
                             Vector{Exonmax}() )

# Main constructor
# Build splice graph here.
function SpliceGraph( gene::Refgene, chrom::DNASequence )
   # splice graph variables
   nodeoffset = Vector{Coordint}()
   nodecoord  = Vector{Coordint}()
   nodelen    = Vector{Coordint}()
   edgetype   = Vector{EdgeType}()
   seq        = DNASequence()
   
   chrom = gene.info[2]
   
   # initialize iterators
   a,alen = 1,length(gene.acc)
   d,dlen = 1,length(gene.don)
   s,slen = 1,length(gene.txst)
   p,plen = 1,length(gene.txen)

   idx = [a, d, s, p]

   curoffset = 1 
  
   while( idx[1] <= alen || idx[2] <= blen || idx[3] <= slen || idx[4] <= plen )    
      # iterate through donors, and acceptors
      # left to right. '-' strand = rc unshift?
      curarr = [ get(gene.txst, idx[1], Inf),
                 get(gene.don,  idx[2], Inf),
                 get(gene.acc,  idx[3], Inf),
                 get(gene.txen, idx[4], Inf) ]

      minidx = indmin( curarr ) 
      minval = min( curarr... )                   
      idx[minidx] += 1

      secarr = [ get(gene.txst, idx[1], Inf),
                 get(gene.don,  idx[2], Inf),
                 get(gene.acc,  idx[3], Inf),
                 get(gene.txen, idx[4], Inf) ]
      secidx = indmin( secarr )
      secval = min( secarr... )

      # last coordinate in the gene:
      if secidx == 1 && getbounds(gene.txst,idx[secidx]) == Inf
         # tack sentinel; break
         break
      end

      # should we make a node?
      if issubinterval( gene.exons, Interval{Coordint}(minval,secval) )
         nodesize = secval - minval #TODO adjustment?
         nodeseq  = chrom[minval:secval] # collect slice
            
      end

      # increment curoffset += nodesize
      
      # determine EdgeType:
      if minidx == 1
         
      elseif minidx == 2

      elseif minidx == 3

      elseif minidx == 4
         
      end

   end # end while

   return SpliceGraph( nodeoffset, nodelen, edgetype, seq )
end

# This function looks specifically for an intersection
# between one interval, and the interval tree, such that
# some interval in the tree full contains the sub interval
#  itree contains: (low,             high)
#  subint must be:     (low+x, high-y)
# Returns: bool
function issubinterval( itree, subint )
   for i in intersect(itree, subint)
      if i.first <= subint.first && i.last >= subint.last
         return true
      end 
   end
   false
end
