
#using Bio.Seq

typealias Str ASCIIString

# overload some stuff here
# lets hash on numeric value for speed 
Base.hash{T <: Bio.Seq.Kmer}( h::T ) = UInt64(h)

# this benchmarks way faster than one initial sizehint! for somereason
function succ_rehash!(d::Dict, targsz::Integer; seed = Int(floor(targsz / 1000)))
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

function bitset!{K}(dict::Dict{K,Bool}, key::K)
   indx = Base.ht_keyindex(dict, key)
   if indx <= 0
     dict[key] = true
   end
end

function __init__()
"""   # This benchmarks faster than BloomFilters.jl
   function set_hash!(d::Dict, n::Integer)
     for i in 1:(n >> 2)
       if i % 1000000 == 0
         println(i)
       end
       k = convert(DNAKmer{32}, UInt64(i))
       d[k] = "Gene1"
     end
   end

   d = Dict{DNAKmer{32},Str}()
   n = 3 * 1000000000 # should be ~ 3billion kmers in the human genome
   @time succ_rehash!(d, n, seed=30000000)
   @time set_hash!(d, n) """
end
