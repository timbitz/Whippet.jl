
using Bio.Seq
using ArgParse

typealias Str ASCIIString

# overload some stuff here
# lets hash on numeric value for speed 
Base.hash{T <: Bio.Seq.Kmer}( h::T ) = UInt64(h)

d = Dict{Kmer,Str}() # set the hash
n = 1000000000 # should be ~ 3billion kmers in the human genome
sizehint!(d, n)

t32 = typeof(DNAKmer(join(["A" for i in 1:32], "")))

# This benchmarks faster than BloomFilters.jl
@time for i in 1:n
  if i % 1000000 == 0
    println(i)
  end
  k = convert(t32, UInt64(i))
  d[k] = "Gene1"
end
