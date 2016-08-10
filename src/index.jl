# index.jl - tim.sterne.weiler@utoronto.ca, 1/28/16

# requires

abstract SeqLibrary

immutable SeqLib <: SeqLibrary
   seq::NucleotideSequence
   offset::Vector{CoordInt}
   names::Vector{GeneName}
   info::Vector{GeneInfo}
   index::FMIndex
   sorted::Bool
end

immutable GraphLib <: SeqLibrary
   offset::Vector{CoordInt}
   names::Vector{GeneName}
   info::Vector{GeneInfo}
   graphs::Vector{SpliceGraph}
   edges::Edges
   index::FMIndex
   sorted::Bool
   kmer::Int64
end

# Binary search.  -- deprecated for more efficient Base.searchsortedlast !!
@inline function search_sorted{T}( arr::Vector{T}, elem::T, low=1, high=length(arr)+1; lower=false )
   low == high && return(lower ? low - 1 : 0)
   const mid = ((high - low) >> 1) + low
   arr[mid] == elem && return(mid)
   if arr[mid] > elem
      return search_sorted(arr, elem, low, mid, lower=lower)
   else
      return search_sorted(arr, elem, mid+1, high, lower=lower)
   end
end

offset_to_name( seqlib::SeqLibrary, offset ) = getindex( seqlib.names, search_sorted(seqlib.offset, CoordInt(offset), lower=true) )

function name_to_offset( seqlib::SeqLibrary, name::GeneName )
   if seqlib.sorted
      ind = search_sorted( seqlib.names, name, true )
      ret = seqlib.names[ind]
   else
      ret = brute_getoffset( seqlib, name )
   end
   ret
end

function brute_getoffset( seqlib::SeqLibrary, name::GeneName )
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

function build_chrom_dict( ref::RefSet )
   ret = Dict{SeqName,Vector{GeneName}}() # refgenomeseq->geneid[]
   for g in keys(ref)
      chrom = ref[g].info.name
      if !haskey(ret, chrom)
         ret[chrom] = Vector{GeneName}()
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
   names = SeqName[]
   for r in fhIter
      immutable!(r.seq)
      push!(names, r.name)
      push!(offset, length(seq)+1)
      println( STDERR, r )
      @time seq *= r.seq
   end
   seq, offset, names
end

function single_genome_index!( fhIter; verbose=false )
   seq,offset,names = load_fasta( fhIter )
   @time fm = FMIndex(twobit_enc(seq), 4, r=4, program=:SuffixArrays, mmap=true)
   println( STDERR, "Finished building Index..." )

   SeqLib(seq, offset, names, fm, false)
end


function trans_index!( fhIter, ref::RefSet; kmer=9 )
   seqdic  = build_chrom_dict( ref )
   xcript  = sg""
   xoffset = Vector{UInt64}()
   xgenes  = Vector{GeneName}()
   xinfo   = Vector{GeneInfo}()
   xgraph  = Vector{SpliceGraph}()
   # set up splice graphs
   runoffset = 0
   for r in fhIter
      haskey(seqdic, r.name) || continue
      Bio.Seq.immutable!(r.seq)
      sg = SGSequence( r.seq )

      r.seq = dna""
      gc() # free

      println(STDERR, "Building Splice Graphs for $( r.name ).." )
      for g in seqdic[r.name]
         #println( STDERR, "Building $g Splice Graph..." )
         curgraph = SpliceGraph( ref[g], sg )
         xcript  *= curgraph.seq
         push!(xgraph, curgraph)
         push!(xgenes, g)
         push!(xinfo, ref[g].info )
         push!(xoffset, runoffset)
         runoffset += length(curgraph.seq) 
      end
   end
   println( STDERR, "Building full sg-index.." )
   @time fm = FMIndex(threebit_enc(xcript), 8, r=1, program=:SuffixArrays, mmap=true) 

   # clean up
   xcript = sg""
   gc()

   println( STDERR, "Building edges.." ) 
   @time edges = build_edges( xgraph, kmer ) # TODO make variable kmer

   GraphLib( xoffset, xgenes, xinfo, xgraph, edges, fm, true, kmer )
end

fixpath( str::String ) = abspath( expanduser( str ) )

function isgzipped( filename::String )
   restr = "\.gz\$"
   re = match(Regex(restr), filename)
   return re == nothing ? false : true
end

function fasta_to_index( filename::String, ref::RefSet; kmer=9 )
   if isgzipped( filename )
      println(STDERR, "Decompressing and Indexing $filename...")
      to_open = open( filename ) |> x->ZlibInflateInputStream(x, reset_on_end=true)
   else
      println(STDERR, "Indexing $filename...")
      to_open = filename
   end
   # iterate through fasta entries
   index = @time trans_index!(open( to_open, FASTA ), ref, kmer=kmer)
   index
end

