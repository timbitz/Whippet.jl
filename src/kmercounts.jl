using Bio.Seq
using FMIndexes
include("index.jl")
@time genome = open(deserialize, "$(pwd())/index/genome.jls")
function kmercounts( genome, k )
   counts = Vector{Int}(4^k)
   for i in 1:(4^k)
      ui = UInt64(i)
      kmer = DNAKmer{k}(ui-1)
      cnt = count(DNASequence(kmer), genome.index)
      counts[i] = cnt
   end
   counts
end
k = 6
c = kmercounts( genome, k )
i = sortperm( c )
DNAKmer{k}(i)
