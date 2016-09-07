function flip( strand::Bio.Intervals.Strand )
    if strand == STRAND_NA || strand == STRAND_BOTH
        return strand
    else
        return convert(Bio.Intervals.Strand, convert(UInt8, strand) $ 0b11 )
    end
end

# "chrX:32049-32077"
function Bio.Intervals.Interval(genome_coord::AbstractString, strand::Union{Bio.Intervals.Strand,Char}=STRAND_BOTH)
    spl_a = split( genome_coord, ':' )
    if length(spl_a) != 2
        error("Invalid genome coord to build Bio.Intervals.Interval!")
    else
        spl_b = split( spl_a[2], '-' )
        if length(spl_b) != 2
           error("Invalid genome coord to build Bio.Intervals.Interval!")
        else
           return Bio.Intervals.Interval(String(spl_a[1]),
                           parse(Int, spl_b[1]),
                           parse(Int, spl_b[2]),
                           strand)
        end
    end
end

function precedes{S, T}(a::Bio.Intervals.Interval{S}, b::Bio.Intervals.Interval{T},
                       seqname_isless::Function=Bio.Intervals.alphanum_isless)
    return (a.last < b.first && a.seqname == b.seqname) ||
        seqname_isless(a.seqname, b.seqname)::Bool
end
