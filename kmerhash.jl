
using Bio.Seq
using ArgParse

typealias Str ASCIIString

# overload some stuff here
# lets hash on identity for speed 
Base.hash{T <: Bio.Seq.Kmer}( h::T ) = UInt64(h)

d = Dict{Kmer,Str}() # set the hash
n = 3000000000 # should be ~ 3billion kmers in the human genome
sizehint!(d, n)

t32 = typeof(DNAKmer(join(["A" for i in 1:32], "")))

# This benchmarks faster than BloomFilters.jl
#@time for i in 1:n
#  k = convert(t32, UInt64(i))
#  d[k] = "Gene1"
#end
