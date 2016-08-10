
typealias SortedOffset Vector{CoordInt}

immutable SeedIndex
   left::Vector{SortedOffset}
   right::Vector{SortedOffset}
   seed_size::Int
end

function SeedIndex( seq, seed_size=18 )
   const kmer_size = seed_size >> 1
   const left  = Vector{SortedOffset}(4^kmer_size)
   const right = Vector{SortedOffset}(4^kmer_size)
   for i in 1:(length(seq) - 2kmer_size + 1)
      lrange,rrange = seedranges( i, kmer_size )
      const lseq = seq[lrange]
      const rseq = seq[rrange]
      try  # this will throw an error if ambiguous nt are present
         const lind = kmer_index(Kmer(lseq))
         const rind = kmer_index(Kmer(rseq))
         isdefinedpush!( left,  lind, convert(CoordInt, i) )
         isdefinedpush!( right, rind, convert(CoordInt, i) )
      end
   end
   sortoffsets!( left )
   sortoffsets!( right )
   SeedIndex( left, right, 2kmer_size )
end 

function seedranges( offset, kmer_size )
   const lrange = offset:(offset+kmer_size-1)
   const rrange = (offset+kmer_size):(offset+2kmer_size-1)
   lrange,rrange
end

@inline function isdefinedpush!( vect::Vector{SortedOffset}, index::Int, value::CoordInt )
   if isdefined( vect, index )
      push!( vect[index], value )
   else
      vect[index] = CoordInt[value]
   end
end

function sortoffsets!( vect::Vector{SortedOffset} )
   for i in 1:length(vect)
      if isdefined( vect, i )
         sort!( vect[i] )
      end
   end
end

@inline function locations( index::SeedIndex, seed )
   if length(seed) != index.seed_size
      error("seed $seed is of incompatible size for SeedIndex!")
   elseif Bio.Seq.find_next_ambiguous( seed, 0 ) > 0
      return CoordInt[]
   end
   const kmer_size = length(seed) >> 1
   const lrange,rrange = seedranges( 1, kmer_size )
   const lind = kmer_index(DNAKmer(seed[lrange]))
   const rind = kmer_index(DNAKmer(seed[rrange]))
   intersect_sorted( index.left[lind], index.right[rind] )
end

@inline function intersect_sorted_null{T}( arrA::Vector{T}, arrB::Vector{T} )
   const res = Nullable{Vector{T}}()
   i,j = 1,1
   while i <= length(arrA) && j <= length(arrB)
      if arrA[i] > arrB[j]
         j += 1
      elseif arrA[i] < arrB[j]
         i += 1
      else
         if isnull(res)
            res = Nullable(Vector{T}(1))
            res.value[1] = arrA[1]
         else
            push!( res.value, arrA[i] )
         end
         i += 1
         j += 1
      end
   end
   res
end
