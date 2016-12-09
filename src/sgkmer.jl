
typealias SGKmer{K} Bio.Seq.Kmer{DNANucleotide, K}

Base.write(io::IO, k::Bio.Seq.Kmer) = write(io, reinterpret(UInt64, k))
Base.read{T <: Bio.Seq.Kmer}(io::IO, ::Type{T}) = convert(T, Base.read(io, UInt64))

# Conversion
# ----------

kmer_index{T,K}( kmer::Bio.Seq.Kmer{T,K} ) = Int(reinterpret(UInt64, kmer)) + 1

kmer_index( seq::BioSequence ) = Int(kmer_index_trailing( seq )) + 1

# calculate kmer index directly
@inline function kmer_index_trailing( seq )
   x     = UInt64(0)
   for nt in seq
      ntint = convert(UInt8, trailing_zeros(nt))
      if ntint > 0x03 || isambiguous(nt)
         return zero(UInt64)
      else
         x = x << 2 | ntint
      end
   end
   x
end

kmer(seq::SGSequence) = sgkmer(seq)
sgkmer(seq::String) = convert(SGKmer, seq)
sgkmer(seq::SGSequence) = convert(SGKmer, seq)
