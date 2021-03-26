# requires:

primitive type EdgeType 8 end

const EDGETYPE_TO_UINT8 = fill( 0x07, (4,4) )
      EDGETYPE_TO_UINT8[4,2] = 0x00 # 'SL' = 0x00; Tx Start
      EDGETYPE_TO_UINT8[4,3] = 0x01 # 'SR' = 0x01; Tandem Last Exon
      EDGETYPE_TO_UINT8[2,4] = 0x02 # 'LS' = 0x02; Tandem First Exon
      EDGETYPE_TO_UINT8[3,4] = 0x03 # 'RS' = 0x03; Tx End
      EDGETYPE_TO_UINT8[2,3] = 0x04 # 'LR' = 0x04; Intron 5'->3' SS
      EDGETYPE_TO_UINT8[2,2] = 0x05 # 'LL' = 0x05; Alt- 5'SS
      EDGETYPE_TO_UINT8[3,3] = 0x06 # 'RR' = 0x06; Alt- 3'SS

const INDEX_TO_EDGETYPE = transpose(reshape([[0x03,0x01,0x03,0x03];
                                             [0x06  for _ in 1:4 ];
                                             [0x05,0x04,0x05,0x02];
                                             [0x00  for _ in 1:4 ] ], (4,4)))

const INDEX_TO_EDGETYPE_NODE = transpose(reshape([[0x03  for _ in 1:4 ];
                                                  [0x06  for _ in 1:4 ];
                                                  [0x05  for _ in 1:4 ];
                                                  [0x00  for _ in 1:4 ] ], (4,4)))


function Base.convert( ::Type{EdgeType}, one::UInt8, two::UInt8 )
   @assert( 5 <= one <= 7 && 5 <= one <= 7 )
   EDGETYPE_TO_UINT8[one-3,two-3]
end

Base.convert( ::Type{EdgeType}, one::DNA, two::DNA ) =
              Base.convert( EdgeType, convert(UInt8, trailing_zeros(one)), convert(UInt8, trailing_zeros(two)) )
Base.convert( ::Type{EdgeType}, edge::UInt8 ) = reinterpret(EdgeType, edge)
Base.convert( ::Type{UInt8}, edge::EdgeType ) = reinterpret(UInt8, edge)
Base.convert( ::Type{I}, edge::EdgeType) where {I <: Integer} = Base.convert(I, Base.convert(UInt8, edge))

EdgeType(edge::UInt8) = reinterpret(EdgeType, edge)
UInt8(edge::EdgeType) = reinterpret(UInt8, edge)

const EDGETYPE_SL = convert(EdgeType, 0b000)
const EDGETYPE_SR = convert(EdgeType, 0b001)
const EDGETYPE_LS = convert(EdgeType, 0b010)
const EDGETYPE_RS = convert(EdgeType, 0b011)

const EDGETYPE_LR = convert(EdgeType, 0b100)
const EDGETYPE_LL = convert(EdgeType, 0b101)
const EDGETYPE_RR = convert(EdgeType, 0b110)

istxstart( edge::EdgeType ) = edge == EDGETYPE_SL || edge == EDGETYPE_LS ? true : false
istxstop(  edge::EdgeType ) = edge == EDGETYPE_SR || edge == EDGETYPE_RS ? true : false

isexonstart( edge::EdgeType ) = edge in (EDGETYPE_LR, EDGETYPE_RR, EDGETYPE_SR, EDGETYPE_LS, EDGETYPE_SL) ? true : false
isexonstop(  edge::EdgeType ) = edge in (EDGETYPE_LR, EDGETYPE_LL, EDGETYPE_SR, EDGETYPE_LS, EDGETYPE_RS) ? true : false
ishardedge(  edge::EdgeType ) = edge in (EDGETYPE_LR, EDGETYPE_LS, EDGETYPE_SR) ? true : false

# Function takes coordinate types for node boundaries
# and returns an EdgeType
function get_edgetype( minidx::Int, secidx::Int, isnode::Bool, strand::Bool )
   if isnode
      ret = INDEX_TO_EDGETYPE_NODE[minidx,secidx]
   else
      ret = INDEX_TO_EDGETYPE[minidx,secidx]
   end
   ret = EdgeType(ret)
   if !strand
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
      return convert(EdgeType, convert(UInt8, edge) ⊻ 0b011)
   end
end

struct GCBiasParams
   length::Int64           # length offset
   width::Float64          # gc bin size
end

const ExpectedGC = Vector{Float16}

function ExpectedGC( seq::BioSequence{A}; gc_param = GCBiasParams(50, 0.05) ) where A <: BioSequences.Alphabet
   bins = zeros(Float16, Int(div(1.0, gc_param.width)+1))
   for i in 1:length(seq)-gc_param.length+1
      gc = Int(div(gc_content(seq[i:i+gc_param.length-1]), gc_param.width)+1)
      @inbounds bins[gc] += 1.0
   end
   # normalize columns
   binsum = sum(bins)
   binsum == 0.0 && (return bins)
   for i in 1:length(bins)
      @fastmath bins[i] /= binsum
   end
   return bins
end

# This holds a representation of the splice graph
# which is a directed multigraph
struct SpliceGraph{K}
   nodeoffset::Vector{CoordInt} # SG offset
   nodecoord::Vector{CoordInt}  # Genome offset
   nodelen::Vector{CoordInt}
   edgetype::Vector{EdgeType}
   edgeleft::Vector{SGKmer{K}}
   edgeright::Vector{SGKmer{K}}
   annopath::Vector{BitSet}
   annoname::Vector{String}
   annobias::Vector{ExpectedGC}
   seq::SGSequence
end
# All positive strand oriented sequences--->
# Node array: txStart| 1 |   2   | 3 |    4    |5| 6 |txEnd
# Edge array:        1   2       3   4         5 6   7
# Node coord:  chr   100 200     300 400      500 600 700

#= empty constructor
SpliceGraph(k::Int) = SpliceGraph( Vector{CoordInt}(), Vector{CoordInt}(),
                                   Vector{CoordInt}(), Vector{EdgeType}(),
                                   Vector{SGKmer{k}}(),Vector{SGKmer{k}}(), dna"" )
=#

# Main constructor
# Build splice graph here.
function SpliceGraph( gene::RefGene, genome::SGSequence, k::Int )
   # splice graph variables
   nodeoffset = Vector{CoordInt}()
   nodecoord  = Vector{CoordInt}()
   nodelen    = Vector{CoordInt}()
   edgetype   = Vector{EdgeType}()
   seq        = LongSequence{DNAAlphabet{4}}()

   strand = gene.info.strand  # Bool

   # initialize iterators
   a,alen = 1,length(gene.acc)
   d,dlen = 1,length(gene.don)
   s,slen = 1,length(gene.txst)
   p,plen = 1,length(gene.txen)

   idx = [a, d, s, p]

   # return a tuple containing the min coordinate's idx,val
   # each of 4 categories, txstart, donor, acceptor, txend
   function getmin_ind_val( gene::RefGene, idx; def=typemax(CoordInt) )
      retarr = CoordInt[ get(gene.txen, idx[1], def),
                         get(gene.acc,  idx[2], def),
                         get(gene.don,  idx[3], def),
                         get(gene.txst, idx[4], def) ]
      argmin( retarr ), min( retarr... )
   end

   while( idx[1] <= alen || idx[2] <= dlen || idx[3] <= slen || idx[4] <= plen )
      # iterate through donors, and acceptors
      # left to right. '-' strand gets revcomp and unshift
      minidx,minval = getmin_ind_val( gene, idx )
      idx[minidx] += 1
      secidx,secval = getmin_ind_val( gene, idx )

      # last coordinate in the gene:
      # if secidx == 1 && get(gene.txst, idx[secidx], Inf) == Inf
      if secval == typemax(CoordInt)
         termedge = strand ? EdgeType(0x03) : EdgeType(0x00)
         stranded_push!(edgetype, termedge, strand)
         break
      end

      # now should we make a node?
      if issubinterval( gene.exons, IntervalTrees.Interval{CoordInt}(minval,secval) ) ||
         (minidx == 2 && minval in gene.novelacc) ||
         (secidx == 3 && secval in gene.noveldon)
         #(usebam && isretainedintron( )) # TODO (not yet implemented)
         leftadj  = (minidx == 1 || minidx == 3) && minval != secval ? 1 : 0
         rightadj = (secidx == 2 || secidx == 4) && minval != secval ? 1 : 0
         nodesize = Int(secval - minval) - leftadj - rightadj + 1
         nodeseq  = copy(genome[(Int(minval)+leftadj):(Int(secval)-rightadj)]) # collect slice
         edge     = get_edgetype( minidx, secidx, true, strand ) # determine EdgeType
         pushval  = minval + leftadj
         thridx = 0
      else # don't make a node, this is a sequence gap, make edge and inc+=2
         idx[secidx] += 1 #skip ahead again
         thridx,thrval = getmin_ind_val( gene, idx )
         rightadj = (thridx == 2 || thridx == 4) && secval != thrval ? 1 : 0
         nodesize = Int(thrval - secval) - rightadj + 1
         nodeseq  = copy(genome[Int(secval):(Int(thrval)-rightadj)])
         edge     = get_edgetype( minidx, secidx, false, strand )
         pushval  = secval
      end

      if strand
         seq *= nodeseq
      else # '-' strand
         seq = reverse_complement(nodeseq) * seq
      end
      stranded_push!(nodecoord, pushval,  strand)
      stranded_push!(nodelen,   nodesize, strand)
      stranded_push!(edgetype,  edge,     strand)

   end # end while

   # calculate node offsets-->
   curoffset = 1
   for n in nodelen
      push!(nodeoffset, curoffset)
      curoffset += n
   end

   eleft  = Vector{SGKmer{k}}(undef, length(edgetype))
   eright = Vector{SGKmer{k}}(undef, length(edgetype))

   paths = build_paths_edges( nodecoord, nodelen, gene )
   names = map( x->x.info.name, gene.reftx )
   expgc = map( x->ExpectedGC(path_to_seq(x, nodeoffset, nodelen, seq)), paths )

   return SpliceGraph( nodeoffset, nodecoord, nodelen, edgetype, eleft, eright, paths, names, expgc, seq )
end

# re-orient - strand by using pushfirst! instead of push!
function stranded_push!( collection, value, strand::Bool )
   if strand
      push!( collection, value )
   else # '-' strand
      pushfirst!( collection, value )
   end
end

# This function looks specifically for an intersection
# between one interval, and the interval tree, such that
# some interval in the tree fully contains the sub interval
#  itree contains: (low,             high)
#  subint must be:     (low+x, high-y)
#  where x and y are >= 0
# Returns: bool
function issubinterval( itree::T, subint::I ) where {T <: IntervalTree,
                                         I <: AbstractInterval}
   for i in intersect(itree, subint)
      if i.first <= subint.first && i.last >= subint.last
         return true
      end
   end
   false
end

# hacky extension of Base.get for tuples which isn't supported natively
# because tuples aren't really a type of collection, fair enough
function Base.get( collection::T, key::I, def ) where {T <: Tuple, I <: Integer}
   if 1 <= key <= length(collection)
      return collection[key]
   else
      return def
   end
end

# Functions for iso.jl annotated edges feature
# Take a RefTx and produce an BitSet through the nodes of a SpliceGraph
# returns: BitSet
function build_annotated_path( nodecoord::Vector{CoordInt},
                               nodelen::Vector{CoordInt},
                               tx::RefTx, strand::Bool )
   path = BitSet()
   # this may be `poor form`, but 256 is too big for default!
   # resize!(path.bits, 64) # Deprecated:  = zeros(UInt32,64>>>5)
   for i in 1:length(tx.acc)
      ind = collect(searchsorted( nodecoord, tx.acc[i], rev=!strand ))[end]
      push!( path, ind )
      cur = ind + (strand ? 1 : -1)
      while 1 <= cur <= length(nodecoord) &&
            nodecoord[cur] <= tx.don[i]
         push!( path, cur )
         cur = cur + (strand ? 1 : -1 )
      end
   end
   path
end

function build_paths_edges( nodecoord::Vector{CoordInt},
                            nodelen::Vector{CoordInt},
                            gene::RefGene )
   paths = Vector{BitSet}()
   for tx in gene.reftx
      curpath = build_annotated_path( nodecoord, nodelen, tx, gene.info.strand )
      push!( paths, curpath )
   end
   paths
end

function path_to_seq( path::BitSet, nodeoffset::Vector{CoordInt},
                      nodelen::Vector{CoordInt}, seq::SGSequence )
   SGSequence(join( map( x->seq[nodeoffset[x]:(nodeoffset[x]+nodelen[x]-1)], collect(path) ) ))
end
