
# This is an overload of the Bio.Seq.complement! function
# to take the complement of all except ambiguous nucleotides
# unlike the original function which complements all nucleotides.
function Bio.Seq.complement!(seq::BioSequence{DNAAlphabet{4}})
    Bio.Seq.orphan!(seq)
    next = Bio.Seq.bitindex(seq, 1)
    stop = Bio.Seq.bitindex(seq, endof(seq) + 1)
    i = 1
    while next < stop
        const ind = Bio.Seq.index(next)
        const x = seq.data[ind]
        seq.data[ind] = (
               ((x & 0x1111111111111111) << 3) | ((x & 0x8888888888888888) >> 3) |
               ((x & 0x2222222222222222) << 1) | ((x & 0x4444444444444444) >> 1))
        if count_ones(x) != 16
           mask = zero(UInt64)
           for j in 0:15
              nt = Bio.Seq.inbounds_getindex(seq, i + j)
              if isambiguous(nt)
                 mask |= 0x000000000000000F << j*4
              end
           end
           seq.data[ind] = (seq.data[ind] & ~mask) | (x & mask)
        end
        next += 64
        i    += 16
    end
    return seq
end

