
using Bio.Seq
using BloomFilters

# overload some stuff here
# lets hash on identity for speed 

Base.hash{T <: Bio.Seq.Kmer}( h::T ) = UInt64(h)

typealias Str ASCIIString

type KmerHash{K <: Kmer, V}
   d::Dict{K,V}
end

bf = BloomFilter(3000000000, 0.01)
