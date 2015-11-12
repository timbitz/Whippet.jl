
using Bio.Seq
using ArgParse


macro nogc(ex)
  quote
    try
      gc_enable(false)
      $(esc(ex))
    finally
      gc_enable(true)
    end
  end
end

typealias Str ASCIIString

# overload some stuff here
# lets hash on numeric value for speed 
Base.hash{T <: Bio.Seq.Kmer}( h::T ) = UInt64(h)

df = Dict{Kmer,Str}() # set the hash
d = Dict{Kmer,Str}()
n = 3000000000 # should be ~ 3billion kmers in the human genome
sizehint!(d, n >> 4)

# This benchmarks faster than BloomFilters.jl
# without sizehint:  16.759637 seconds (263.34 M allocations: 6.512 GB, 27.96% gc time)
# with sizehint:
function set_hash!(d::Dict, n::Integer)
  for i in 1:(n >> 5)
    if i % 1000000 == 0
      println(i)
    end
    
    k = convert(DNAKmer{32}, UInt64(i))
    d[k] = "Gene1"
  end
  nothing
end

@time @nogc set_hash!(d, n)
@time @nogc set_hash!(df, n)
