
typealias SGKmer{K} Bio.Seq.Kmer{DNANucleotide, K}

Base.write(io::IO, k::Bio.Seq.Kmer) = write(io, reinterpret(UInt64, k))
Base.read{T <: Bio.Seq.Kmer}(io::IO, ::Type{T}) = convert(T, Base.read(io, UInt64))

# Conversion
# ----------

kmer_index{T,K}( kmer::Bio.Seq.Kmer{T,K} ) = Int(UInt64(kmer)) + 1

kmer_index( seq::BioSequence ) = Int(kmer_index_trailing( seq )) + 1
#kmer_index( seq::SGSequence  ) = Int(kmer_index_straight( seq )) + 1

# calculate kmer index directly
function kmer_index_trailing( seq )
   x     = UInt64(0)
   for nt in seq
      ntint = convert(UInt8, trailing_zeros(nt))
      if ntint > 0x03 | isambiguous(nt)
         return 0
      else
         x = x << 2 | ntint
      end
   end
   x
end

# calculate kmer index directly
#=function kmer_index_straight( seq )
   x     = UInt64(0)
   for nt in seq
      ntint = convert(UInt8, nt)
      if ntint > 0x03 | isambiguous(nt)
         return 0
      else
         x = x << 2 | ntint
      end
   end
   x
end

Base.convert(::Type{SGKmer}, seq::SGSequence) = sgkmer(seq)

kmer(seq::SGSequence) = sgkmer(seq)
sgkmer(seq::String) = convert(SGKmer, seq)
function sgkmer(seq::SGSequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(SGKmer{length(seq)}, kmer_index_straight(seq))
end=#
