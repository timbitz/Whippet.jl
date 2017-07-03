@inline function raw_getindex{A}(seq::BioSequence{A}, i::Integer)
    j = Bio.Seq.bitindex(seq, i)
    @inbounds return convert(UInt8, (seq.data[Bio.Seq.index(j)] >> Bio.Seq.offset(j)) & Bio.Seq.mask(A))
end


## This function is necessary as a result of changes made in Bio.jl as of Julia v0.5
## Since all nucleotides are converted to UInt8 as one hot 4-bit encodings by default
## we can override FMIndexes.sa_range function for those sequences specifically
## to re-encode as DNAAlphabet{2}
@inline function FMIndexes.sa_range{A <: BioSequences.Alphabet}(query::BioSequence{A}, index::FMIndex, init_range::UnitRange{Int})
    sp, ep = init_range.start, init_range.stop
    # backward search
    i = length(query)
    @inbounds while sp ≤ ep && i ≥ 1
        char = raw_getindex(query, i)
        c = index.count[char+1]
        sp = c + rank(char, index.bwt, (sp > index.sentinel ? sp - 1 : sp) - 1) + 1
        ep = c + rank(char, index.bwt, (ep > index.sentinel ? ep - 1 : ep))
        i -= 1
    end
    return sp:ep
end
