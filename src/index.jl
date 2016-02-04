# index.jl - tim.sterne.weiler@utoronto.ca, 1/28/16

# requires
using FMIndexes
using IntArrays

include("types.jl")
include("bio_nuc_safepatch.jl")
include("refflat.jl")
include("graph.jl")
include("edges.jl")

typealias Str ASCIIString

abstract SeqLibrary

immutable SeqLib <: SeqLibrary
   seq::NucleotideSequence
   offset::Vector{Coordint}
   names::Vector{Str}
   index::FMIndex
   sorted::Bool
end

immutable GraphLib <: SeqLibrary
   offset::Vector{Coordint}
   names::Vector{Str}
   graphs::Vector{SpliceGraph}
   #edges::Vector{Edges}
   index::FMIndex
   sorted::Bool
end

function getoffset( seqlib::SeqLibrary, name::Str )
   if seqlib.sorted
      ret = sorted_getoffset( seqlib, name )
   else
      ret = brute_getoffset( seqlib, name )
   end
   ret
end

function sorted_getoffset( seqlib, name, low=1, high=length(seqlib.names)+1 )
   low == high && return(-1)
   mid = ((high - low) >> 1) + lowe
   seqlib.names[mid] == name && return(seqlib.offset[mid])
   if seqlib.names[mid] > name
      ret = sorted_getoffset(seqlib, name, low, mid)
   else
      ret = sorted_getoffset(seqlib, name, mid+1, high)
   end
   ret
end

function brute_getoffset( seqlib::SeqLibrary, name::Str )
   ret = -1
   for n in 1:length(seqlib.names)
     if seqlib.names[n] == name
       ret = seqlib.offset[n]
       break
     end
   end
   ret
end

function build_offset_dict{I <: Integer, 
                           S <: AbstractString}( offset::Vector{I}, names::Vector{S} )
   ret = Dict{S,I}()
   for i in 1:length(names)
      cname,coffset = names[i],offset[i]
      ret[cname] = coffset
   end
   ret
end

function build_chrom_dict( ref::Refset )
   ret = Dict{Str,Vector{Genename}}() # refgenomeseq->geneid[]
   for g in keys(ref.geneset)
      chrom = ref.geneset[g].info[1]
      if !haskey(ret, chrom)
         ret[chrom] = Vector{Genename}()
      end
      push!(ret[chrom], g)
   end
   ret
end

# replace L,R,S,N with A
function twobit_enc(seq)
   len = length(seq)
   ret = IntVector{2,UInt8}(len)
   for i in 1:len
      if seq[i] == DNA_N
         ret[i] = convert(UInt8, DNA_A)
      else
         ret[i] = convert(UInt8, seq[i])
      end
   end
   ret
end

# encode 3-bit sequence with L,R,S,N,A,T,G,C
function threebit_enc(seq)
   len = length(seq)
   ret = IntVector{3,UInt8}(len)
   for i in 1:len
      ret[i] = convert(UInt8, seq[i])
   end
   ret
end

function load_fasta( fhIter; verbose=false )
   seq = SGSequence(mutable=false)
   offset = Int[]
   names = Str[]
   for r in fhIter
      immutable!(r.seq)
      push!(names, r.name)
      push!(offset, length(seq)+1)
      println( STDERR, r )
      @time seq *= r.seq
   end
   seq, offset, names
end

#TODO test offset
function load_fasta_sg( fhIter; verbose=false )
   seq = sg""
   offset = Int[]
   names = Str[]
   for r in fhIter
      Bio.Seq.immutable!(r.seq)
      sg = SGSequence( r.seq )
      push!(names, r.name)
      push!(offset, length(seq)+1)
      println( STDERR, r )
      @time seq *= sg 
   end    
   seq, offset, names
end

function single_genome_index!( fhIter; verbose=false )
   seq,offset,names = load_fasta( fhIter )
   @time fm = FMIndex(twobit_enc(seq), 4, r=4, program=:SuffixArrays, mmap=true)
   println( STDERR, "Finished building Index..." )

   SeqLib(seq, offset, names, fm, false)
end


function trans_index!( fhIter, ref::Refset )
   seqdic  = build_chrom_dict( ref )
   xcript  = sg""
   xoffset = Vector{UInt64}()
   xgenes  = Vector{Genename}()
   xgraph  = Vector{SpliceGraph}()
   # set up splice graphs
   runoffset = 0
   for r in fhIter
      Bio.Seq.immutable!(r.seq)
      sg = SGSequence( r.seq )

      r.seq = dna""
      gc() # free

      haskey(seqdic, r.name) || continue
      for g in seqdic[r.name]
         println( STDERR, "Building $g Splice Graph..." )
         curgraph = SpliceGraph( ref.geneset[g], sg )
         xcript  *= curgraph.seq
         push!(xgraph, curgraph)
         push!(xgenes, g)
         push!(xoffset, runoffset)
         runoffset += length(curgraph.seq) 
      end
   end
   @time fm = FMIndex(threebit_enc(xcript), 8, r=1, program=:SuffixArrays, mmap=true) 
   println( STDERR, "Finished building index..." )

   # clean up
   xcript = sg""
   gc()

   GraphLib( xoffset, xgenes, xgraph, fm, true)
end

