using IntervalTrees


typealias Str ASCIIString

immutable Seqlibrary
   seq::DNASequence
   offset::Vector{Int64}
   names::Vector{Str}
   index::FMIndex
end

function getoffset( genome::Genome, name::Str )
   ret = -1
   for n in 1:length(genome.names)
     if genome.names[n] == name
       ret = genome.offset[n]
       break
     end
   end
   ret
end

# replace N with A
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


function threebit_enc(seq)
   len = length(seq)
   ret = IntVector{3,UInt8}(len)
   for i in 1:len
      ret[i] = convert(UInt8, seq[i])
   end
   ret
end

function load_fasta( fhIter; verbose=false )
   seq = DNASequence()
   offset = Int[]
   names = Str[]
   immutable!(seq)
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
   Seqlibrary(seq, offset, names, fm)
end

function trans_index!( fhIter, ref::Refflat )
   seq,offset,names = load_fasta( fhIter )
   for g in keys(ref.genetotx)
      
   end
   Seqlibrary( )
end

