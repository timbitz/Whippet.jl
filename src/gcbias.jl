function bitparallel_gc_content( seq::BioSequence, start=1, stop=length(seq) )
    gc   = 0
    bs = bitindex(seq, start)
    be = bitindex(seq, stop)
    for i in index(bs):index(be)
        @inbounds x = seq.data[i]
        # mask nucleotides upstream of first(seq.part)
        if i == index(bs)
            x &= 0xFFFFFFFFFFFFFFFF << offset(bs)
        end
        # mask nucleotides downstream of last(seq.part)
        if i == index(be)
            x &= 0xFFFFFFFFFFFFFFFF >> (64 - (offset(be)+4))
        end
        # count GC nucleotides
        a =  x & 0x1111111111111111
        c = (x & 0x2222222222222222) >> 1
        g = (x & 0x4444444444444444) >> 2
        t = (x & 0x8888888888888888) >> 3
        gc += count_ones((c | g) & ~(a | t))
    end
    return length(seq) > 0 ? gc / (stop-start+1) : 0.0
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
            gc = div(bitparallel_gc_content(seq, i, i+len-1), gcp.width)
            @inbounds data[gc, col] += 1.0
        end
    end
    # normalize columns
    for c in 1:size(bins, 2)
        bins[:,c] /= sum(bins[:,c])
    end
    bins
end


