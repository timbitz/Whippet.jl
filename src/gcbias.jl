
function parallel_gc_content( seq::BioSequence, start=1, stop=length(seq) )
    gc   = 0
    s = bitindex(seq, 1)
    e = bitindex(seq, length(seq))
    for i in index(s):index(e)
        @inbounds data = seq.data[i]
        # mask nucleotides upstream of first(seq.part)
        if i == index(s)
            data &= 0xFFFFFFFFFFFFFFFF << offset(s)
        end
        # mask nucleotides downstream of last(seq.part)
        if i == index(e)
            data &= 0xFFFFFFFFFFFFFFFF >> (64 - (offset(e)+4))
        end
        # count unambiguous GC nucleotides
        gc += count_ones(data & 0x6666666666666666)
    end
    return length(seq) > 0 ? gc / length(seq) : 0.0
end

struct GCBiasParams
    offset::Int64           # length offset
    increment::Int64        # length increment
    max::Int64              # length maximum
    width::Float64          # gc bin size
end

const ExpectedGC = Array{Float16,2}

function ExpectedGC( seq::BioSequence; gc_param = GCBiasParams(offset=50, 
                                                               increment=25, 
                                                               max=convert(Int,typemax(UInt8)), 
                                                               width=0.05) )

    bins = zeros(Float16, div(1.0, gcp.width), div(gcp.max-gcp.offset, gcp.increment))
    col  = 1
    for len in gcp.offset:gcp.increment:gcp.max
        for i in 1:length(seq)-len
            gc = div(parallel_gc_content(seq, i, i+len-1), gcp.width)
            @inbounds data[gc, col] += 1.0
        end
    end
    # normalize columns
    for c in 1:size(bins, 2)
        bins[:,c] /= sum(bins[:,c])
    end
    bins
end


