#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Bio.Seq

include("kmerhash.jl")

K = 32
N = 3 * 100000000 # should be ~ 3billion kmers in the human genome

println(STDERR, "Allocating k-mer hash...")
gdict = Dict{DNAKmer{K},Bool}()
@time succ_rehash!(gdict, N, seed=3000000)

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
         println(STDERR, r )
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
open(string("gdict", ".index"), "w+") do io
    serialize( io, gdict )
end


