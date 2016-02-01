# requires:
# include("bio_nuc_patch.jl")
using IntervalTrees

bitstype 8 EdgeType

const EDGETYPE_TO_UINT8 = fill( 0x07, (4,4) )
EDGETYPE_TO_UINT8[4,2] = 0x00 # 'SL' = 0x00; Tx Start 
EDGETYPE_TO_UINT8[4,3] = 0x01 # 'SR' = 0x01; Tandem Last Exon
EDGETYPE_TO_UINT8[2,4] = 0x02 # 'LS' = 0x02; Tandem First Exon
EDGETYPE_TO_UINT8[3,4] = 0x03 # 'RS' = 0x03; Tx End
EDGETYPE_TO_UINT8[2,3] = 0x04 # 'LR' = 0x04; Intron 5'->3' SS
EDGETYPE_TO_UINT8[2,2] = 0x05 # 'LL' = 0x05; Alt- 5'SS
EDGETYPE_TO_UINT8[3,3] = 0x06 # 'RR' = 0x06; Alt- 3'SS

const INDEX_TO_EDGETYPE = transpose(reshape([[0x00  for _ in 1:4 ]; 
                                             [0x02,0x05,0x04,0x05];
                                             [0x06  for _ in 1:4 ];
                                             [0x03,0x03,0x01,0x03] ], (4,4)))

const INDEX_TO_EDGETYPE_NODE = transpose(reshape([[0x00  for _ in 1:4 ];
                                                  [0x05  for _ in 1:4 ];
                                                  [0x06  for _ in 1:4 ];
                                                  [0x03  for _ in 1:4 ] ], (4,4)))

const EDGETYPE_TO_DNA = DNASequence[ dna"SL", dna"SR", dna"LS", dna"RS",
                                     dna"LR", dna"LL", dna"RR", dna"SS" ]

function Base.convert( ::Type{EdgeType}, one::UInt8, two::UInt8 )
   @assert( 5 <= one <= 7 && 5 <= one <= 7 ) 
   EDGETYPE_TO_UINT8[one-3,two-3]
end

Base.convert( ::Type{EdgeType}, one::DNANucleotide, two::DNANucleotide ) = 
               Base.convert( EdgeType, convert(UInt8, one), convert(UInt8, two) )
Base.convert( ::Type{EdgeType}, edge::UInt8 ) = box(EdgeType, unbox(UInt8, edge ))
Base.convert( ::Type{UInt8}, edge::EdgeType ) = box(UInt8, unbox(EdgeType, edge ))
Base.convert{I <: Integer}( ::Type{I}, edge::EdgeType) = Base.convert(I, Base.convert(UInt8, edge))
Base.convert( ::Type{DNASequence}, edge::EdgeType ) = EDGETYPE_TO_DNA[Base.convert(UInt8, edge)+1]

# Function takes coordinate types for node boundaries
# and returns an EdgeType
function get_edgetype( minidx::Int, secidx::Int, isnode::Bool, strand::Char='+' )
   if isnode
      ret = INDEX_TO_EDGETYPE_NODE[minidx,secidx]
   else
      ret = INDEX_TO_EDGETYPE[minidx,secidx]
   end
   ret = EdgeType(ret)
   if strand != '+'
      ret = invert_edgetype( ret )
   end
   ret
end

# Invert EdgeTypes using xor
# 'LR' does not need to be inverted
function invert_edgetype( edge::EdgeType )
   if edge == EdgeType(0x04) # 'LR'
      return edge
   else
      return convert(EdgeType, convert(UInt8, edge) $ 0b011)
   end
end

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
   seq        = dna""
   
   strand = gene.info[2]
   
   # initialize iterators
   a,alen = 1,length(gene.acc)
   d,dlen = 1,length(gene.don)
   s,slen = 1,length(gene.txst)
   p,plen = 1,length(gene.txen)

   idx = [a, d, s, p]

   # return a tuple containing the min coordinate's idx,val
   # each of 4 categories, txstart, donor, acceptor, txend
   function getmin_ind_val( gene::Refgene, idx; def=Inf )
      retarr = [ get(gene.txst, idx[1], def),
                 get(gene.don,  idx[2], def),
                 get(gene.acc,  idx[3], def),
                 get(gene.txen, idx[4], def) ]
      indmin( retarr ), min( retarr... )
   end

   # re-orient - strand by using unshift! instead of push!
   function stranded_push!( collection, value, strand::Char )
      if strand == '+'
         push!( collection, value )
      else # '-' strand
         unshift!( collection, value )
      end
   end

   while( idx[1] <= alen || idx[2] <= dlen || idx[3] <= slen || idx[4] <= plen )    
      # iterate through donors, and acceptors
      # left to right. '-' strand = rc unshift?
      minidx,minval = getmin_ind_val( gene, idx )
      idx[minidx] += 1
      secidx,secval = getmin_ind_val( gene, idx )

      # last coordinate in the gene:
#      if secidx == 1 && get(gene.txst, idx[secidx], Inf) == Inf
      if secval == Inf
         termedge = EdgeType(0x03)
         stranded_push!(edgetype, termedge, strand)
         if strand == '+'
            seq *= DNASequence(termedge)
         else
            termedge = invert_edgetype( termedge )
            seq = DNASequence(termedge) * seq
         end
         break
      end
      
      # now should we make a node?
      if issubinterval( gene.exons, Interval{Coordint}(minval,secval) )
         nodesize = secval - minval #TODO adjustment?
         nodeseq  = dna"AAAAA" #chrom[minval:secval] # collect slice
         edge = get_edgetype( minidx, secidx, true, strand ) # determine EdgeType
         pushval = minval
         thridx = 0
      else # don't make a node, this is a sequence gap, make edge and inc+=2
         idx[secidx] += 1 #skip ahead again
         thridx,thrval = getmin_ind_val( gene, idx )
         nodesize = thrval - secval
         nodeseq = dna"CCCCC"
         edge = get_edgetype( minidx, secidx, false, strand )
         pushval = secval
      end

      if strand == '+'
         seq *= DNASequence(edge) * nodeseq
      else # '-' strand
         seq = reverse_complement(nodeseq) * DNASequence(edge) * seq
      end
      println("strand: $strand, minidx: $minidx, secidx: $secidx, thridx: $thridx, nodesize: $nodesize, edgetype: $edge, pushval: $(Int(pushval)), nodeseq: $nodeseq")
      stranded_push!(nodecoord, pushval,  strand)
      stranded_push!(nodelen,   nodesize, strand)
      stranded_push!(edgetype,  edge,     strand)

   end # end while
   
   # need to calculate nodeoffsets now.

   return SpliceGraph( nodeoffset, nodecoord, nodelen, edgetype, seq )
end

# This function looks specifically for an intersection
# between one interval, and the interval tree, such that
# some interval in the tree full contains the sub interval
#  itree contains: (low,             high)
#  subint must be:     (low+x, high-y)
#  where x and y are >= 0
# Returns: bool
function issubinterval{T <: IntervalTree, 
                       I <: AbstractInterval}( itree::T, subint::I )
   for i in intersect(itree, subint)
      if i.first <= subint.first && i.last >= subint.last
         return true
      end 
   end
   false
end

# hacky extension of Base.get for tuples which isn't supported natively
# because tuples aren't really a type of collection, fair enough
function Base.get{T <: Tuple, I <: Integer}( collection::T, key::I, def )
   if 1 <= key <= length(collection)
      return collection[key]
   else
      return def
   end
end
