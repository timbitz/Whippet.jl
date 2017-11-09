
const SGKmer{K} = BioSequences.Kmer{DNA, K}

Base.write(io::IO, k::BioSequences.Kmer) = write(io, reinterpret(UInt64, k))
Base.read{T <: BioSequences.Kmer}(io::IO, ::Type{T}) = convert(T, Base.read(io, UInt64))

# Conversion
# ----------

kmer_index{T,K}( kmer::BioSequences.Kmer{T,K} ) = Int(reinterpret(UInt64, kmer)) + 1

kmer_index( seq::BioSequence ) = Int(kmer_index_trailing( seq )) + 1

# calculate kmer index directly
@inline function kmer_index_trailing( ::Type{T}, seq ) where T <: Unsigned
   x     = T(0)
   for nt in seq
      ntint = convert(UInt8, trailing_zeros(nt))
      if ntint > 0x03 || isambiguous(nt)
         return zero(T)
      else
         x = x << 2 | ntint
      end
   end
   x
end

kmer_index_trailing( seq ) = kmer_index_trailing(UInt64, seq )

kmer(seq::SGSequence) = sgkmer(seq)
sgkmer(seq::String) = convert(SGKmer, seq)
sgkmer(seq::SGSequence) = convert(SGKmer, seq)
