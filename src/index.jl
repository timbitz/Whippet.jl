using IntervalTrees

typealias Str ASCIIString

immutable Seqlibrary
   seq::DNASequence
   offset::Vector{Int64}
   names::Vector{Str}
   index::FMIndex
   sorted::Bool
end

function getoffset( seqlib::Seqlibrary, name::Str )
   if seqlib.sorted
      ret = sorted_getoffset( seqlib, name )
   else
      ret = brute_getoffset( seqlib, name )
   end
   ret
end

function sorted_getoffset( seqlib, name, low=1, high=length(seqlib.names)+1 )
   low == high && return(-1)
   mid = ((high - low) >> 1) + low
   seqlib.names[mid] == name && return(seqlib.offset[mid])
   if seqlib.names[mid] > name
      ret = sorted_getoffset(seqlib, name, low, mid)
   else
      ret = sorted_getoffset(seqlib, name, mid+1, high)
   end
   ret
end

function brute_getoffset( seqlib::Seqlibrary, name::Str )
   ret = -1
   for n in 1:length(seqlib.names)
     if seqlib.names[n] == name
       ret = seqlib.offset[n]
       break
     end
   end
   ret
end

# replace J, S, N with A
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

# encode 3-bit sequence with J,S,N,A,T,G,C
function threebit_enc(seq)
   len = length(seq)
   ret = IntVector{3,UInt8}(len)
   for i in 1:len
      ret[i] = convert(UInt8, seq[i])
   end
   ret
end

function load_fasta( fhIter; verbose=false )
   seq = DNASequence(mutable=false)
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

function single_genome_index!( fhIter; verbose=false )
   seq,offset,names = load_fasta( fhIter )
   @time fm = FMIndex(twobit_enc(seq), 4, r=4, program=:SuffixArrays, mmap=true)
   println( STDERR, "Finished building Index..." )
   Seqlibrary(seq, offset, names, fm, false)
end

function trans_index!( fhIter, ref::Refset )
   seq,offset,names = load_fasta( fhIter )
   xcript  = DNASequence(mutable=false)
   xoffset = Vector{UInt64}()
   xgenes  = Vector{Genename}()
   genes = sort( collect(keys(ref.gninfo)), by= k->ref.gninfo[k][1] ) # by chrom
   for g in genes
      # set up exons
      # insert L/R and S, push!         
   end
   #Seqlibrary( , true)
end

