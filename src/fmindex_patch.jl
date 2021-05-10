using WaveletMatrices

for n in (2, 4)
    @eval begin
        bitsof(::Type{DNAAlphabet{$n}}) = $n
        bitsof(::Type{RNAAlphabet{$n}}) = $n
    end
end

@inline function raw_getindex(seq::BioSequence{A}, i::Integer) where A
    j = BioSequences.bitindex(seq, i)
    @inbounds return convert(UInt8, (seq.data[BioSequences.index(j)] >> BioSequences.offset(j)) & BioSequences.bitmask(bitsof(A)))
end

## This function is necessary as a result of changes made in Bio.jl as of Julia v0.5
## Since all nucleotides are converted to UInt8 as one hot 4-bit encodings by default
## we can override FMIndexes.sa_range function for those sequences specifically
## to re-encode as DNAAlphabet{2}
@inline function FMIndexes.sa_range(query::BioSequence{A}, index::FMIndex, init_range::UnitRange{Int}) where A <: BioSequences.Alphabet
    sp, ep = init_range.start, init_range.stop
    # backward search
    i = length(query)
    @inbounds while sp ≤ ep && i ≥ 1
        char = raw_getindex(query, i)
        c = index.count[char+1]
        sp = c + WaveletMatrices.rank(char, index.bwt, (sp > index.sentinel ? sp - 1 : sp) - 1) + 1
        ep = c + WaveletMatrices.rank(char, index.bwt, (ep > index.sentinel ? ep - 1 : ep))
        i -= 1
    end
    return sp:ep
end
