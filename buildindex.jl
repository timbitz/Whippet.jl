#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Bio.Seq
using FMIndexes
using IntArrays
using GZip

include("kmerhash.jl")

K = 32
N = 100000000 # should be ~ 3billion kmers in the human genome

println(STDERR, "Allocating k-mer hash...")
#gdict = Dict{DNAKmer{K},Bool}()
#@time succ_rehash!(gdict, N, seed=10000000)
gfmind = FMIndex[]

# encode a DNA sequence with 3-bit unsigned integers;
# this is because a reference genome has five nucleotides: A/C/G/T/N.
function encode(seq)
    encoded = IntVector{3,UInt8}(length(seq))
    for i in 1:endof(seq)
        encoded[i] = convert(UInt8, seq[i])
    end
    return encoded
end

for f in readdir()
   re = match(r"(\S+).(fa|fasta)(.gz)?", f)
  
   if re != nothing
      mainname = re.captures[1]
      if re.captures[3] != nothing #then gzipped
         println(STDERR, "Decompressing and Indexing $mainname...")
         #to_open = pipeline(`cat $f`, `gzip -dc`)
         continue 
      else
         println(STDERR, "Indexing $mainname...")
         to_open = f
      end
      # iterate through fasta entries
      fh = open( to_open, FASTA )
      for r in fh
         immutable!(r.seq)
         println(STDERR, r )
         @time fm = FMIndex(encode(r.seq))
         push!(gfmind, fm)
         r.name == "chr2" ? break : continue
         @time for (i,k) in each(DNAKmer{K}, r.seq, K >> 1)
            bitset!(gdict, k)
            if i % 262145 == 0
               print(STDERR, "Hashing Kmer $i...\r")
            end
         end
      end
   end 
end
println()
println(STDERR, "Saving genome index")
GZip.open(string("gfmind", ".jls"), "w+") do io
   @time serialize( io, gfmind )
end


