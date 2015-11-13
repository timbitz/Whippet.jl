
using Bio.Seq
using ArgParse

typealias Str ASCIIString

# overload some stuff here
# lets hash on numeric value for speed 
Base.hash{T <: Bio.Seq.Kmer}( h::T ) = UInt64(h)

df = Dict{DNAKmer{32},Str}() # set the hash
d = Dict{DNAKmer{32},Str}()
n = 3000000000 # should be ~ 3billion kmers in the human genome
@time succ_rehash!(d, n)

function succ_rehash!(d::Dict, targsz::Integer)
  seed = 1000000
  if targsz <= seed
    return
  else
    cursz = seed
    while cursz <= targsz
      Base.rehash!(d, cursz)
      cursz *= 10
    end
  end
end

# This benchmarks faster than BloomFilters.jl
# without sizehint:  16.759637 seconds (263.34 M allocations: 6.512 GB, 27.96% gc time)
# with sizehint:
function set_hash!(d::Dict, n::Integer)
  for i in 1:(n >> 4)
    if i % 1000000 == 0
      println(i)
    end
    k = convert(DNAKmer{32}, UInt64(i))
    d[k] = "Gene1"
  end
end

@time set_hash!(d, n)
@time set_hash!(df, n)
