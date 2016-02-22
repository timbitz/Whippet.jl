#  This is a port of the source code from https://github.com/BioJulia/Bio.jl/blob/master/src/seq/nucleotide.jl
#  This file is written almost entirely by the authors of BioJulia/Bio.jl
#  It has been altered here to include new types containing masked nucleotides 
#  'L' for donor upstream splice site
#  'R' for acceptor downstream splice site, and 
#  'S' for transcription start/stop site.  
#  These can then be encoded with 'N' and 'A','T','C','G' into 3-bits

#@everywhere using Bio
#using Bio.Seq

using Base.Intrinsics
import Base.length,Base.start,Base.*,Base.^,Base.done,Base.==,Base.next,Base.reverse

#abstract Sequence
# Nucleotides
# ===========


# Single nucleotides are represented in bytes using just the two low-order bits
#abstract Nucleotide
bitstype 8 SGNucleotide <: Bio.Seq.Nucleotide

# Conversion from/to integers
# ---------------------------

Base.convert(::Type{SGNucleotide}, nt::UInt8) = box(SGNucleotide, unbox(UInt8, nt))
Base.convert(::Type{UInt8}, nt::SGNucleotide) = box(UInt8, unbox(SGNucleotide, nt))
#Base.convert{T<:Unsigned, S<:Nucleotide}(::Type{T}, nt::S) = box(T, Base.zext_int(T, unbox(S, nt)))
#Base.convert{T<:Unsigned, S<:Nucleotide}(::Type{S}, nt::T) = convert(S, convert(UInt8, nt))


# Nucleotide encoding definition
# ------------------------------

# SG Nucleotides

"SG Adenine"
const SG_A = convert(SGNucleotide, 0b000)

"SG Cytosine"
const SG_C = convert(SGNucleotide, 0b001)

"SG Guanine"
const SG_G = convert(SGNucleotide, 0b010)

"SG Thymine"
const SG_T = convert(SGNucleotide, 0b011)

"SG Any Nucleotide"
const SG_N = convert(SGNucleotide, 0b100)

"SG Invalid Nucleotide"
const SG_INVALID = convert(SGNucleotide, 0b1000) # Indicates invalid SG when converting string

## Patch
"SG Left Edge Signal"
const SG_L = convert(SGNucleotide, 0b101)

"SG Right Edge Signal"
const SG_R = convert(SGNucleotide, 0b110)

"SG Start/Stop Signal"
const SG_S = convert(SGNucleotide, 0b111) 

"Returns Any SG Nucleotide (SG_N)"
nnucleotide(::Type{SGNucleotide}) = SG_N
invalid_nucleotide(::Type{SGNucleotide}) = SG_INVALID

## Patch
lnucleotide(::Type{SGNucleotide}) = SG_L
rnucleotide(::Type{SGNucleotide}) = SG_R
snucleotide(::Type{SGNucleotide}) = SG_S

lnucleotide{N <: Nucleotide}(::Type{N}) = SG_L
rnucleotide{N <: Nucleotide}(::Type{N}) = SG_R
snucleotide{N <: Nucleotide}(::Type{N}) = SG_S

function isvalid(nt::SGNucleotide)
    return convert(UInt8, nt) ≤ convert(UInt8, SG_S)
end



# Conversion from Char
# --------------------

# lookup table for characters in 'A':'t'
const char_to_sg = [
    SG_A,       SG_INVALID, SG_C,       SG_INVALID, SG_INVALID, SG_INVALID,
    SG_G,       SG_INVALID, SG_INVALID, SG_INVALID, SG_INVALID, SG_L,
    SG_INVALID, SG_N,       SG_INVALID, SG_INVALID, SG_INVALID, SG_R,
    SG_S,        SG_T,       SG_INVALID, SG_INVALID, SG_INVALID, SG_INVALID,
    SG_INVALID, SG_INVALID, SG_INVALID, SG_INVALID, SG_INVALID, SG_INVALID,
    SG_INVALID, SG_INVALID, SG_A,       SG_INVALID, SG_C,       SG_INVALID,
    SG_INVALID, SG_INVALID, SG_G,       SG_INVALID, SG_INVALID, SG_INVALID,
    SG_INVALID, SG_L,        SG_INVALID, SG_N,       SG_INVALID, SG_INVALID,
    SG_INVALID, SG_R,        SG_S,        SG_T ]

@inline function Base.convert(::Type{SGNucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 't' ? char_to_sg[c - 'A' + 1] : SG_INVALID
    return nt
end

@inline function unsafe_ascii_byte_to_nucleotide(T::Type{SGNucleotide}, c::UInt8)
    @inbounds nt = char_to_sg[c - 0x41 + 1]
    return nt
end



# Conversion to Char
# ------------------

##
const sg_to_char = ['A', 'C', 'G', 'T', 'N', 'L', 'R', 'S']

Base.convert(::Type{Char}, nt::SGNucleotide) = sg_to_char[convert(UInt8, nt) + 1]

# Converstion to DNANucleotide
# ----------------------------
==( a::Bio.Seq.DNANucleotide, b::SGNucleotide ) = box(UInt8, a) == box(UInt8, b)
==( a::SGNucleotide, b::Bio.Seq.DNANucleotide ) = box(UInt8, a) == box(UInt8, b)


# Basic functions
# ---------------

function Base.show(io::IO, nt::SGNucleotide)
    if !isvalid(nt)
        write(io, "Invalid SG Nucleotide")
    else
        write(io, convert(Char, nt))
    end
end



# Nucleotide Sequence
# ===================

# A general purpose nucleotide sequence representation
#
# Nucleotide sequences are 2-bit encoded and packed into UInt64s. 'N's are
# represented with an N mask BitArray. If the ns[i] bit is set, then the
# sequence may have any nucleotide at that position and it must be ignored.
#
# Sequences are either explicitly immutable or mutable. Immutable sequences have
# the benefit that subsequence operations (`seq[i:j]`) are cheap, since they
# share the underlying data and do not copy anything. Mutable sequences may be
# mutated, but as a consequence subsequences copy data rather than reference it.
#
# Sequences can be converted between mutable and immutable using `mutable!` and
# `immutable!`, respectively. Converting from mutable to immutable is cheap: it
# only flips the `mutable` flag. The converse, immutable to mutable, is cheap
# *if* the sequence is not a subsequence and has no subsequences, otherwise it
# must make a copy of the data.
#
##
type NucleotideSequence{T<:SGNucleotide} <: Bio.Seq.Sequence
    data::Vector{UInt64} # 2-bit encoded sequence
    ns::BitVector        # 'N' mask
    ls::BitVector        # 'L' donor splice site
    rs::BitVector        # 'R' acceptor splice site
    ss::BitVector        # 'S' mask
    part::UnitRange{Int} # interval within `data` and `ns` defining the (sub)sequence
    mutable::Bool        # true if the sequence can be safely mutated

    # true if this was constructed as a subsequence of another sequence or if
    # subsequences were constructed from this sequence. When this is true, we
    # need to copy the data to convert from immutable to mutable
    hasrelatives::Bool
end

# Constructors
# ------------

"""
`NucleotideSequence(SGNucleotide|DEPRECATEDNucleotide)`

Construct an empty nucleotide sequence of the given type
"""
NucleotideSequence{T<:SGNucleotide}(::Type{T}; mutable::Bool=false) =
    NucleotideSequence{T}(zeros(UInt64, 0), BitVector(0),  BitVector(0),  BitVector(0), 1:0, mutable, false)


"""
`NucleotideSequence(SGNucleotide|DEPRECATEDNucleotide, other::NucleotideSequence, part::UnitRange)`

Construct a subsequence of the given type from another nucleotide sequence
"""
function NucleotideSequence{T<:SGNucleotide}(::Type{T}, other::NucleotideSequence,
                                           part::UnitRange; mutable::Bool=false)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    if start < other.part.start || stop > other.part.stop
        error("Invalid subsequence range")
    end

    seq = NucleotideSequence{T}(other.data, other.ns, other.ls, other.rs, other.ss, start:stop, mutable, true)

    if other.mutable || mutable
        orphan!(seq, true)
    end

    if !other.mutable
        other.hasrelatives = true
    end

    if !mutable
        seq.hasrelatives = true
    end

    return seq
end


# Faster divrem(n, 32)
function divrem32(n::Integer)
    return (n >> 5, n & 0b11111)
end


# Faster divrem(n, 64)
function divrem64(n::Integer)
    return (n >> 6, n & 0b111111)
end


# Return the number of UInt64s needed to represent a sequence of length n
function seq_data_len(n::Integer)
    return cld(n, 32)
end

# This is the body of the NucleotideSequence constructor below. It's separated
# into a macro so we can generate two versions depending on wether the `unsafe`
# flag is set.
macro encode_seq(nt_convert_expr, strdata, seqdata, ns, ls, rs, ss)
    quote
        j = startpos
        idx = 1
        @inbounds begin
            # we OR all the nucleotides to detect and if 0b1000 is set,
            # then we know there was an invalid nucleotide.
            ored_nucs = UInt8(0)

            for i in 1:length($(seqdata))
                shift = 0
                data_i = UInt64(0)
                while shift < 64 && j <= stoppos
                    c = $(strdata)[j]
                    nt = $(nt_convert_expr)
                    if nt == nnucleotide(T)
                        # manually inlined: ns[i] = true
                        d = (idx - 1) >>> 6
                        r = (idx - 1) & 63
                        $(ns).chunks[d + 1] |= UInt64(1) << r
                    elseif nt == lnucleotide(T)
                        # manually inlined: ns[i] = true
                        d = (idx - 1) >>> 6
                        r = (idx - 1) & 63
                        $(ls).chunks[d + 1] |= UInt64(1) << r
                    elseif nt == rnucleotide(T)
                        # manually inlined: ns[i] = true
                        d = (idx - 1) >>> 6
                        r = (idx - 1) & 63
                        $(rs).chunks[d + 1] |= UInt64(1) << r
                    elseif nt == snucleotide(T)
                        # manually inlined: ns[i] = true
                        d = (idx - 1) >>> 6
                        r = (idx - 1) & 63
                        $(ss).chunks[d + 1] |= UInt64(1) << r
                    else
                        ored_nucs |= convert(UInt8, nt)
                        data_i |= convert(UInt64, nt) << shift
                    end

                    j += 1
                    idx += 1
                    shift += 2
                end
                $(seqdata)[i] = data_i
            end

            if ored_nucs & 0b1000 != 0
                # invalid nucleotide: figure out what the first bad character was.
                for i in startpos:stoppos
                    c = $(strdata)[j]
                    nt = $(nt_convert_expr)
                    if nt == invalid_nucleotide(T)
                        error(string(c, " is not a valid ", T))
                    end
                end
            end
        end
    end
end


"""
`NucleotideSequence(SGNucleotide|DEPRECATEDNucleotide, seq::AbstractString, startpos::Int, stoppos::Int)`

Construct a nucleotide sequence from the `seq[startpos:stoppos]` string
"""
function NucleotideSequence{T<:SGNucleotide}(::Type{T}, seq::Union{AbstractString, Vector{UInt8}},
                                           startpos::Int, stoppos::Int,
                                           unsafe::Bool=false; mutable::Bool=false)
    len = stoppos - startpos + 1
    data = zeros(UInt64, seq_data_len(len))
    ns = falses(len)
    rs = falses(len)
    ls = falses(len)
    ss = falses(len)

    if unsafe
        @encode_seq(unsafe_ascii_byte_to_nucleotide(T, c), seq, data, ns, ls, rs, ss)
    else
        @encode_seq(convert(T, convert(Char, c)), seq, data, ns, ls, rs, ss)
    end

    return NucleotideSequence{T}(data, ns, ls, rs, ss, 1:len, mutable, false)
end

function NucleotideSequence{T<:SGNucleotide}(t::Type{T}, seq::Union{AbstractString, Vector{UInt8}}; mutable::Bool=false)
    return NucleotideSequence(t, seq, 1, length(seq), mutable=mutable)
end

"""
`NucleotideSequence(seq::AbstractVector{T<:SGNucleotide}, startpos::Int, stoppos::Int)`

Construct a nucleotide sequence from the `seq[startpos:stoppos]` vector
"""
function NucleotideSequence{T<:SGNucleotide}(seq::AbstractVector{T},
                                           startpos::Int, stoppos::Int,
                                           unsafe::Bool=false; mutable::Bool=false)
    len = stoppos - startpos + 1
    data = zeros(UInt64, seq_data_len(len))
    ns = falses(len)
    ls = falses(len)
    rs = falses(len)
    ss = falses(len)

    if unsafe
        @encode_seq(c, seq, data, ns, ls, rs, ss)
    else
        @encode_seq(begin
            if !isvalid(c)
                error("the sequence includes a valid nucleotide at $j")
            end
            c
        end, seq, data, ns, ls, rs, ss)
    end

    return NucleotideSequence{T}(data, ns, ls, rs, ss, 1:len, mutable, false)
end


"""
Reset the contents of a mutable sequence from a string.
"""
function Base.copy!{T<:SGNucleotide}(seq::NucleotideSequence{T}, strdata::Vector{UInt8},
                  startpos::Integer, stoppos::Integer)
    if !seq.mutable
        error("Cannot copy! to immutable sequnce. Call `mutable!(seq)` first.")
    end

    n = stoppos - startpos + 1
    len = seq_data_len(n)
    if length(seq.data) < len
        resize!(seq.data, len)
    end

    if length(seq.ns) < n
        resize!(seq.ns, n)
    end

    if length(seq.ls) < n
        resize!(seq.ls, n)
    end

    if length(seq.rs) < n
        resize!(seq.rs, n)
    end

    if length(seq.ss) < n
        resize!(seq.ss, n)
    end

    fill!(seq.data, 0)
    fill!(seq.ns, false)
    fill!(seq.ls, false)
    fill!(seq.rs, false)
    fill!(seq.ss, false)
    seq.part = 1:n

    @encode_seq(convert(T, convert(Char, c)), strdata, seq.data, seq.ns, seq.ls, seq.rs, seq.ss)

    return seq
end


"""
`NucleotideSequence(chunks::NucleotideSequence...)`

Construct a nucleotide sequence by concatenating the given sequences
"""
function NucleotideSequence{T<:SGNucleotide}(chunks::NucleotideSequence{T}...)
    seqlen = 0
    for chunk in chunks
        seqlen += length(chunk)
    end

    datalen = seq_data_len(seqlen)
    data = zeros(UInt64, datalen)
    ns   = falses(seqlen)
    ls   = falses(seqlen)
    rs   = falses(seqlen)
    ss   = falses(seqlen)
    newseq = NucleotideSequence{T}(data, ns, ls, rs, ss, 1:seqlen, false, false)

    pos = 1
    for chunk in chunks
        unsafe_copy!(newseq, pos, chunk)
        pos += length(chunk)
    end

    return newseq
end


(*){T}(chunk1::NucleotideSequence{T}, chunks::NucleotideSequence{T}...) = NucleotideSequence(chunk1, chunks...)


"""
Construct a NucleotideSequence from an array of nucleotides.
"""
function NucleotideSequence{T<:SGNucleotide}(seq::AbstractVector{T}; mutable::Bool=false)
    return NucleotideSequence(seq, 1, endof(seq), mutable=mutable)
end


# Mutability/Immutability
# -----------------------

ismutable(seq::NucleotideSequence) = seq.mutable


function mutable!(seq::NucleotideSequence)
    if !seq.mutable
        if seq.hasrelatives
            orphan!(seq, true)
        end
        seq.mutable = true
    end
    return seq
end


function immutable!(seq::NucleotideSequence)
    seq.mutable = false
    return seq
end


"""
`repeat(chunk::NucleotideSequence, n)`

Construct a nucleotide sequence by repeating another sequence `n` times
"""
function repeat{T<:SGNucleotide}(chunk::NucleotideSequence{T}, n::Integer)
    seqlen = n * length(chunk)

    datalen = seq_data_len(seqlen)
    data = zeros(UInt64, datalen)
    ns   = falses(seqlen)
    ls   = falses(seqlen)
    rs   = falses(seqlen)
    ss   = falses(seqlen)
    newseq = NucleotideSequence{T}(data, ns, ls, rs, ss, 1:seqlen, false, false)

    pos = 1
    for i in 1:n
        unsafe_copy!(newseq, pos, chunk)
        pos += length(chunk)
    end

    return newseq
end

"Repeat nucleotide sequences"
(^){T}(chunk::NucleotideSequence{T}, n::Integer) = repeat(chunk, n::Integer)


"""
Copy `src` to `dest` starting at position `pos`.

This is unsafe in the following ways:
- Disregards immutability of `dest`
- May write a few bases past `dest[pos + length(src) - 1]`
- Doesn't bounds check anything.

It's really only suitable for use in the concatenation constructor.
"""
## TODO
function Base.unsafe_copy!{T}(dest::NucleotideSequence{T}, pos::Int, src::NucleotideSequence{T})
    abspos = dest.part.start + pos - 1
    copy!(dest.ns, abspos, src.ns, src.part.start, length(src))
    copy!(dest.ls, abspos, src.ls, src.part.start, length(src))
    copy!(dest.rs, abspos, src.rs, src.part.start, length(src))
    copy!(dest.ss, abspos, src.ss, src.part.start, length(src))

    d1, r1 = divrem(abspos - 1, 32)
    d2, r2 = divrem(src.part.start - 1, 32)

    l = 0
    while l < length(src)
        if r1 == r2 == 0
            dest.data[d1+1] = src.data[d2+1]
        else
            dest.data[d1+1] |= (src.data[d2+1] >> (2*r2)) << (2*r1)
        end

        if r1 > r2
            k = 32 - r1
            d1 += 1
            r1 = 0
            r2 += k
            l += k
        elseif r2 > r1
            k = 32 - r2
            d2 += 1
            r2 = 0
            r1 += k
            l += k
        else # r1 == r2
            k = 32 - r1
            r1 += k
            r2 += k
            if r1 >= 32
                r1 = 0
                r2 = 0
                d1 += 1
                d2 += 1
            end
            l += k
        end
    end

    # zero positions that we've overwritten at the end
    if l > length(src)
        if r1 == 0
            r1 = 32
            d1 -= 1
        end
        dest.data[d1+1] &= 0xffffffffffffffff >> (2 * (32 - r1 + l - length(src)))
    end
end


# Aliases and contructors
# -----------------------

# SG Sequences
typealias SGSequence NucleotideSequence{SGNucleotide}

"Construct an empty SG nucleotide sequence"
SGSequence(; mutable::Bool=true) =
    NucleotideSequence(SGNucleotide, mutable=mutable)

"Construct a SG nucleotide subsequence from another sequence"
SGSequence(other::NucleotideSequence, part::UnitRange; mutable::Bool=false) =
    NucleotideSequence(SGNucleotide, other, part, mutable=mutable)

"Construct a SG nucleotide sequence from an AbstractString"
SGSequence(seq::AbstractString; mutable=false) =
    NucleotideSequence(SGNucleotide, seq, mutable=mutable)

"Construct a SG nucleotide sequence from other sequences"
SGSequence(chunk1::SGSequence, chunks::SGSequence...) = NucleotideSequence(chunk1, chunks...)
SGSequence(seq::Union{Vector{UInt8}, AbstractString}; mutable::Bool=false) =
    NucleotideSequence(SGNucleotide, seq, mutable=mutable)
SGSequence(seq::Union{Vector{UInt8}, AbstractString}, startpos::Int, endpos::Int, unsafe::Bool=false; mutable::Bool=false) =
    NucleotideSequence(SGNucleotide, seq, startpos, endpos, unsafe, mutable=mutable)
SGSequence(seq::AbstractVector{SGNucleotide}; mutable::Bool=false) =
    NucleotideSequence(seq, mutable=mutable)



# Conversion
# ----------

# Convert from/to vectors of nucleotides
Base.convert{T<:SGNucleotide}(::Type{NucleotideSequence},    seq::AbstractVector{T}) = NucleotideSequence(seq, 1, endof(seq))
Base.convert{T<:SGNucleotide}(::Type{NucleotideSequence{T}}, seq::AbstractVector{T}) = NucleotideSequence(seq, 1, endof(seq))
Base.convert{T<:SGNucleotide}(::Type{Vector{T}}, seq::NucleotideSequence{T}) = [x for x in seq]

# Convert from/to Strings
Base.convert(::Type{SGSequence}, seq::AbstractString) = SGSequence(seq)
Base.convert(::Type{AbstractString}, seq::NucleotideSequence) = convert(ASCIIString, [convert(Char, x) for x in seq])

# Convert between DNA/RNA and SG
Base.convert(::Type{SGSequence}, seq::Bio.Seq.NucleotideSequence) = SGSequence( convert(AbstractString, seq) )



# Basic Functions
# ---------------

Base.length(seq::NucleotideSequence) = length(seq.part)
Base.endof(seq::NucleotideSequence)  = length(seq)

# Pretting printing of sequences.
function Base.show{T}(io::IO, seq::NucleotideSequence{T})
    len = length(seq)
    write(io, "$(string(len))nt ",
          seq.mutable ? "Mutable " : "",
          T == SGNucleotide ? "SG" : "DEPRECATED", " Sequence\n ")

    # don't show more than this many characters to avoid filling the screen
    # with junk
    const maxcount = Base.tty_size()[2] - 2
    if len > maxcount
        for nt in seq[1:div(maxcount, 2) - 1]
            write(io, convert(Char, nt))
        end
        write(io, "…")
        for nt in seq[(end - (div(maxcount, 2) - 1)):end]
            write(io, convert(Char, nt))
        end
    else
        for nt in seq
            write(io, convert(Char, nt))
        end
    end
end

function =={T}(a::NucleotideSequence{T}, b::NucleotideSequence{T})
    if a.data === b.data && a.ns === b.ns && a.part == b.part
        return true
    end

    if length(a) != length(b)
        return false
    end

    for (u, v) in zip(a, b)
        if u != v
            return false
        end
    end

    return true
end

# Get the nucleotide at position i, ignoring the N mask.
@inline function getnuc(T::Type, data::Vector{UInt64}, i::Integer)
    d, r = divrem32(i - 1)
    return convert(T, convert(UInt8, (data[d + 1] >>> (2*r)) & 0b11))
end

# This function returns true if there are no n's in the sequence
# it returns false if there are. TODO use parts. for now use hasn?
function iskmersafe{T <: Bio.Seq.NucleotideSequence}( seq::T )
   for i in 1:length( seq.ns.chunks )
      seq.ns.chunks[i] == 0 || return false
   end
   return true
end

function Base.getindex{T}(seq::NucleotideSequence{T}, i::Integer)
    if i > length(seq) || i < 1
        throw(BoundsError())
    end
    i += seq.part.start - 1
    if seq.ns[i]
        return nnucleotide(T)
    elseif seq.ls[i]
        return lnucleotide(T)
    elseif seq.rs[i]
        return rnucleotide(T)
    elseif seq.ss[i]
        return snucleotide(T)
    else
        return getnuc(T, seq.data, i)
    end
end

# Construct a subesequence
Base.getindex{T}(seq::NucleotideSequence{T}, r::UnitRange) = NucleotideSequence{T}(seq, r)

Base.setindex!{T}(seq::NucleotideSequence{T}, r::UnitRange) = NucleotideSequence{T}(seq, r)


function Base.setindex!{T<:SGNucleotide}(seq::NucleotideSequence{T}, nt::T, i::Integer)
    if !seq.mutable
        error("Cannot mutate an immutable sequence. Call `mutable!(seq)` first.")
    end
    if i < 1 || i > length(seq)
        throw(BoundsError())
    end

    d, r = divrem(i - 1, 32)
    if nt == nnucleotide(T)
        @inbounds seq.ns[i] = true
        @inbounds seq.data[d + 1] $= UInt64(0b11) << (2*r)
    elseif nt == lnucleotide(T)
        @inbounds seq.ls[i] = true
        @inbounds seq.data[d + 1] $= UInt64(0b11) << (2*r)
    elseif nt == rnucleotide(T)
        @inbounds seq.rs[i] = true
        @inbounds seq.data[d + 1] $= UInt64(0b11) << (2*r)
    elseif nt == snucleotide(T)
        @inbounds seq.ss[i] = true
        @inbounds seq.data[d + 1] $= UInt64(0b11) << (2*r)
    else
        @inbounds seq.ns[i] = false
        @inbounds seq.ls[i] = false
        @inbounds seq.rs[i] = false
        @inbounds seq.ss[i] = false 
        @inbounds seq.data[d + 1] =
            (seq.data[d + 1] & ~(UInt64(0b11) << (2*r))) |
                (convert(UInt64, nt) << (2*r))
    end
    return nt
end

Base.setindex!{T}(seq::NucleotideSequence{T}, nt::Char, i::Integer) =
    setindex!(seq, convert(T, nt), i)


# Replace a NucleotideSequence's data with a copy, copying only what's needed.
#
# If a sequence's starting position within its data != 1, this function
# will copy a subset of the data to align with sequences range and make start == 1.
#
# The user should never need to call this, as it has no outward effect on the
# sequence, but it makes functions like mismatch easier and faster if can assume
# a sequence is aligned with its data.
#
# If reorphan = true, copy the data regardless of the start position.
#
function orphan!{T<:SGNucleotide}(seq::NucleotideSequence{T}, reorphan=false)
    if !reorphan && seq.part.start == 1
        return seq
    end

    data   = zeros(UInt64, seq_data_len(length(seq)))
    d0, r0 = divrem32(seq.part.start - 1)

    h = 64 - 2*r0
    k = 2*r0

    j = d0 + 1
    @inbounds for i in 1:length(data)
        data[i] |= seq.data[j] >>> k

        j += 1
        if j > length(seq.data)
            break
        end

        data[i] |= seq.data[j] << h
    end

    seq.data = data
    seq.ns   = seq.ns[seq.part.start:seq.part.stop]
    seq.ls   = seq.ls[seq.part.start:seq.part.stop]
    seq.rs   = seq.rs[seq.part.start:seq.part.stop]
    seq.ss   = seq.ss[seq.part.start:seq.part.stop]
    seq.part = 1:length(seq.part)
    seq.hasrelatives = false

    return seq
end

# Copy a sequence.
#
# Unlike constructing subsequences with seq[a:b], this function actually copies
# the underlying data. Since sequences are immutable, you should basically
# never have to do this. It's useful only for low-level algorithms like
# `reverse` which actually do make a copy and modify the copy in
# place.
Base.copy{T<:SGNucleotide}(seq::NucleotideSequence{T}) =
    orphan!(NucleotideSequence{T}(seq.data, seq.ns, seq.ls, seq.rs, seq.ss, seq.part, seq.mutable, seq.hasrelatives), true)


# Iterating throug nucleotide sequences
@inline function Base.start{T<:SGNucleotide}(seq::NucleotideSequence{T})
    npos = nextone(seq.ns, seq.part.start)
    lpos = nextone(seq.ls, seq.part.start)
    rpos = nextone(seq.rs, seq.part.start)
    spos = nextone(seq.ss, seq.part.start)
    return seq.part.start, npos, lpos, rpos, spos
end

@inline function Base.next{T<:SGNucleotide}(seq::NucleotideSequence{T}, state)
    i, npos,lpos,rpos,spos = state
    if i == npos
        npos = nextone(seq.ns, i + 1)
        value = nnucleotide(T)
    elseif i == lpos
        lpos = nextone(seq.ls, i + 1)
        value = lnucleotide(T)
    elseif i == rpos
        rpos = nextone(seq.rs, i + 1)
        value = rnucleotide(T)
    elseif i == spos
        spos = nextone(seq.ss, i + 1)
        value = snucleotide(T)
    else
        d, r = divrem32(i - 1)
        @inbounds value = convert(T, ((seq.data[d + 1] >>> (2 * r)) & 0b11) % UInt8)
    end

    return value, (i + 1, npos, lpos, rpos, spos)
end

@inline Base.done(seq::NucleotideSequence, state) = state[1] > seq.part.stop

# String Decorator
# ----------------

# Enable building sequence literals like: sg"ACGTACGT" and deprecated"ACGUACGU"
macro sg_str(seq, flags...)
    return SGSequence(remove_newlines(seq))
end

function remove_newlines(seq)
    return replace(seq, r"\r|\n", "")
end


# Transformations
# ---------------

# In-place complement (nucleotide complement is equivalent to bitwise complement
# in the encoding used)
function unsafe_complement!{T<:SGNucleotide}(seq::NucleotideSequence{T})
    @inbounds for i in 1:length(seq.data)
        seq.data[i] = ~seq.data[i]
    end
    return seq
end


"""
`complement(seq::NucleotideSequence)`

The nucleotide complement of the sequence `seq`
"""
Base.complement{T<:SGNucleotide}(seq::NucleotideSequence{T}) = unsafe_complement!(copy(seq))

# Nucleotide reverse. Reverse a kmer stored in a UInt64.
function nucrev(x::UInt64)
     x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>>  2
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >>> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >>> 32
     return x
end

"""
`reverse(seq::NucleotideSequence)`

Reversed copy of the nucleotide sequence `seq`
"""
function Base.reverse{T<:SGNucleotide}(seq::NucleotideSequence{T})
    if isempty(seq)
        return seq
    end

    orphan!(seq)

    k = (2 * length(seq) + 63) % 64 + 1
    h = 64 - k

    data = zeros(UInt64, length(seq.data))
    j = length(data)
    i = 1
    @inbounds while true
        x = nucrev(seq.data[i])
        data[j] |= x >>> h
        if (j -= 1) == 0
            break
        end
        data[j] |= x << k;
        i += 1
    end

    return NucleotideSequence{T}(data, reverse(seq.ns), reverse(seq.ls), reverse(seq.rs), reverse(seq.ss), seq.part, seq.mutable, false)
end

# Return the reverse complement of seq
# Return a reversed copy of seq
"""
`reverse_complement(seq::NucleotideSequence)`

Reversed complement of the nucleotide sequence `seq`
"""
reverse_complement{T<:SGNucleotide}(seq::NucleotideSequence{T}) = unsafe_complement!(reverse(copy(seq)))


# SequenceNIterator
# -----------------

# Iterate through positions in the sequence with Ns
#
# This can be much faster than testing every position (seq.ns[i]) since
# it can skip over 64 positions at time if they don't have 'N's.
immutable SequenceNIterator
    ns::BitVector
    part::UnitRange{Int}
end

SequenceNIterator(seq::NucleotideSequence) = SequenceNIterator(seq.ns, seq.part)
npositions(seq::NucleotideSequence)        = SequenceNIterator(seq)


# Find the next N in the sequence starting at position i.
#
# Return any position past the end of the sequence if there are no more Ns.
#
function nextn(it::SequenceNIterator, i)
    d, r = divrem64(i - 1)
    while d < length(it.ns.chunks) && it.ns.chunks[d + 1] >>> r == 0 && d * 64 < it.part.stop
        d += 1
        r  = 0
    end

    if d * 64 + r + 1 > it.part.stop
        return d * 64 + r + 1
    end

    if d + 1 <= length(it.ns.chunks)
        x = it.ns.chunks[d + 1] >>> r
        while x & 0x1 == 0
            x >>>= 1
            r += 1
        end
    end

    return d * 64 + r + 1
end


# Find the index of the next 1 bit, starting at index i.
function nextone(b::BitVector, i)
    if i > length(b)
        return i
    end

    d, r = divrem64(i - 1)
    chunk = b.chunks[d + 1] >>> r
    if chunk != 0
        t = trailing_zeros(chunk)
        return i + t
    end

    i += 64 - r
    r = 0
    d += 1
    @inbounds for l in (d+1):length(b.chunks)
        if b.chunks[l] != 0
            return i + trailing_zeros(b.chunks[l])
        end
        i += 64
    end

    return min(length(b) + 1, i)
end


# Iterating through SequenceNIterator
# -----------------------------------

Base.start(it::SequenceNIterator) = nextn(it, it.part.start)

function Base.next(it::SequenceNIterator, i)
    d, r = divrem64(i - 1)
    next_i = nextn(it, i + 1)
    return i + it.part.start - 1, next_i
end

Base.done(it::SequenceNIterator, i) = i > it.part.stop

function hasn(seq::NucleotideSequence)
    it = npositions(seq)
    return !done(it, start(it))
end

# TODO: Implement length for SequenceNIterators to use comprehensions


# Mismatch counting
# -----------------

"Mismatch count between two kmers"
function nucmismatches(x::UInt64, y::UInt64)
    xyxor = x $ y
    return count_ones((xyxor & 0x5555555555555555) | ((xyxor & 0xAAAAAAAAAAAAAAAA) >>> 1))
end

"Mask of the first `k` bits of a UInt64"
function makemask(k::Integer)
    return 0xffffffffffffffff >> (64 - k)
end

"""
`mismatches(a::NucleotideSequence, b::NucleotideSequence, [nmatches=false])`

Return the number of mismatches between `a` and `b`.

If `a` and `b` are of differing lengths, only the first `min(length(a), length(b))`
nucleotides are compared.

### Arguments
  * `a`: first sequence to compare
  * `b`: second sequence to compare
  * `nmatches`: if true, N matches anything, if false, N matches only itself (false)

### Returns
The number of mismatches
"""
#=
function mismatches{T}(a::NucleotideSequence{T}, b::NucleotideSequence{T},
                       nmatches::Bool=false)

    # we need to assume that a is aligned with its data, rearrange the
    # comparison if that's not the case.
    if a.part.start != 1
        if b.part.start == 1
            return mismatches(b, a)
        else
            if length(a) < b
                orphan!(a)
            else
                orphan!(b)
                return mismatches(b, a)
            end
        end
    end

    # Count mismatches, ignoring the presence of 'N's in the sequence.
    d0, r0 = divrem64(b.part.start - 1)
    h = 64 - r0
    k = r0

    hmask = makemask(h)
    kmask = makemask(k)
    count = 0
    j = d0 + 1
    for i in 1:length(a.data)
        count += nucmismatches(a.data[i] & hmask, b.data[j] >>> k)
        j += 1
        if j > length(b.data)
            break
        end
        count += nucmismatches(a.data[i] & kmask, b.data[j] << h)
    end

    # Here's the ugly part. We've just counted mismtaches without taking the N
    # mask into account. If 'N's are present in the sequence, that mismatch
    # count may be too high or too low, so we walk through all N positions in
    # both sequences in unision, correcting the count where needed.
    nsa       = npositions(a)
    nsb       = npositions(b)
    nsa_state = start(nsa)
    nsb_state = start(nsb)
    a_done    = done(nsa, nsa_state)
    b_done    = done(nsb, nsb_state)

    i, nsa_state = a_done ? (0, nsa_state) : next(nsa, nsa_state)
    j, nsb_state = b_done ? (0, nsb_state) : next(nsb, nsb_state)

    while true
        a_hasnext = !done(nsa, nsa_state)
        b_hasnext = !done(nsb, nsb_state)

        if !a_done && (b_done || i < j)
            nucsmatch = getnuc(T, a.data, i + a.part.start - 1) ==
                        getnuc(T, b.data, i + b.part.start - 1)
            if nucsmatch
                if !nmatches
                    count += 1
                end
            else
                if nmatches
                    count -= 1
                end
            end

            if a_hasnext; i, nsa_state = next(nsa, nsa_state)
            else; a_done = true; end
        elseif !b_done && (a_done || j < i)
            nucsmatch = getnuc(T, a.data, j + a.part.start - 1) ==
                        getnuc(T, b.data, j + b.part.start - 1)
            if nucsmatch
                if !nmatches
                    count += 1
                end
            else
                if nmatches
                    count -= 1
                end
            end

            if b_hasnext; j, nsb_state = next(nsb, nsb_state)
            else; b_done = true; end
        elseif !a_done && !b_done
            nucsmatch = getnuc(T, a.data, i + a.part.start - 1) ==
                        getnuc(T, b.data, i + b.part.start - 1)
            if !nucsmatch
                count -= 1
            end

            if a_hasnext; i, nsa_state = next(nsa, nsa_state)
            else; a_done = true; end
            if b_hasnext; j, nsb_state = next(nsb, nsb_state)
            else; b_done = true; end
        else
            break
        end
    end

    return count
end
=#

# Additional read/write functions
Base.write(io::Base.TCPSocket, nt::Bio.Seq.QualityEncoding) = Base.write(io, convert(UInt16, nt))
Base.read{T <: Bio.Seq.QualityEncoding}(io::Base.TCPSocket, t::Type{T}) = convert(T, Base.read(io, UInt16))

# K-mer
# =====


# A Kmer is a sequence <= 32nt, without any 'N's, packed in a single 64 bit value.
#
# While NucleotideSequence is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in DEPRECATED-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.

bitstype 64 Kmer{T<:Nucleotide, K}

typealias SGKmer{K} Kmer{SGNucleotide, K}



# Conversion
# ----------

# Conversion to/from UInt64

Base.convert{K}(::Type{SGKmer{K}}, x::UInt64) = box(SGKmer{K}, unbox(UInt64, x))
Base.convert{K}(::Type{UInt64}, x::SGKmer{K}) = box(UInt64, unbox(SGKmer{K}, x))


Base.write{T <: Bio.Seq.Kmer}(io::Base.IOStream, k::T) = Base.write(io, convert(UInt64, k))
Base.write{T <: Kmer}(io::Base.IOStream, k::T) = Base.write(io, convert(UInt64, k))
Base.write{T <: Kmer}(io::Base.TCPSocket, k::T) = Base.write(io, convert(UInt64, k))

Base.read{T <: Bio.Seq.Kmer}(io::Base.IOStream, t::Type{T}) = convert(T, Base.read(io, UInt64))
Base.read{T <: Kmer}(io::Base.IOStream, t::Type{T}) = convert(T, Base.read(io, UInt64))
Base.read{T <: Kmer}(io::Base.TCPSocket, t::Type{T}) = convert(T, Base.read(io, UInt64))

==( a::Bio.Seq.DNAKmer, b::SGKmer ) = convert(UInt64, a) == convert(UInt64, b)
==( a::SGKmer, b::Bio.Seq.DNAKmer ) = convert(UInt64, a) == convert(UInt64, b)

# Conversion to/from String

function Base.convert{T, K}(::Type{Kmer{T, K}}, seq::AbstractString)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    @assert length(seq) == K error("Cannot construct a $(K)-mer from a string of length $(length(seq))")

    x     = UInt64(0)
    shift = 0
    for (i, c) in enumerate(seq)
        nt = convert(T, c)
        @assert nt != nnucleotide(T) error("A K-mer may not contain an N in its sequence")
        @assert nt != lnucleotide(T) error("A K-mer may not contain an L in its sequence")
        @assert nt != rnucleotide(T) error("A K-mer may not contain an R in its sequence")
        @assert nt != snucleotide(T) error("A K-mer may not contain an S in its sequence")

        x |= convert(UInt64, nt) << shift
        shift += 2
    end

    return convert(Kmer{T, K}, x)
end

Base.convert{T}(::Type{Kmer{T}}, seq::AbstractString) = convert(Kmer{T, length(seq)}, seq)
Base.convert{T, K}(::Type{AbstractString}, seq::Kmer{T, K}) = convert(AbstractString, [convert(Char, x) for x in seq])

# Conversion to/from NucleotideSequence

"Convert a NucleotideSequence to an SGKmer"
function Base.convert{T, K}(::Type{Kmer{T,K}}, seq::NucleotideSequence{T})
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    @assert length(seq) == K error("Cannot construct a $(K)-mer from a NucleotideSequence of length $(length(seq))")

    x     = UInt64(0)
    shift = 0
    for (i, nt) in enumerate(seq)
        if nt == nnucleotide(T) || nt == lnucleotide(T) || nt == rnucleotide(T) || nt == snucleotide(T)
            error("A K-mer may not contain an N,L,R,S in its sequence")
        end
        x |= convert(UInt64, nt) << shift
        shift += 2
    end
    return convert(Kmer{T, K}, x)
end

Base.convert{T}(::Type{Kmer},    seq::NucleotideSequence{T}) = convert(Kmer{T, length(seq)}, seq)
Base.convert{T}(::Type{Kmer{T}}, seq::NucleotideSequence{T}) = convert(Kmer{T, length(seq)}, seq)

function Base.convert{T, K}(::Type{NucleotideSequence{T}}, x::Kmer{T, K})
    ns = falses(K)
    return NucleotideSequence{T}([convert(UInt64, x)], ns, 1:K, false, false)
end
Base.convert{T, K}(::Type{NucleotideSequence}, x::Kmer{T, K}) = convert(NucleotideSequence{T}, x)


==( a::Bio.Seq.DNANucleotide, b::SGNucleotide ) = convert(UInt64, a) == convert(UInt64, b)

# Constructors
# ------------

# From strings

"Construct a SGKmer to an AbstractString"
sgkmer(seq::AbstractString) = convert(SGKmer, seq)

"Construct a Kmer from a sequence of Nucleotides"
function kmer{T <: Nucleotide}(nts::T...)
    K = length(nts)
    if K > 32
        error(string("Cannot construct a K-mer longer than 32nt."))
    end

    x = UInt64(0)
    shift = 0
    for (i, nt) in enumerate(nts)
        if nt == nnucleotide(T)
            error("A Kmer may not contain an N in its sequence")
        end
        x |= convert(UInt64, nt) << shift
        shift += 2
    end
    return convert(Kmer{T, K}, x)
end

"Construct a Kmer from a SGSequence"
function kmer(seq::SGSequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(SGKmer{length(seq)}, seq)
end


# call kmer with @inline macro would reduce the performance significantly?
# Would the compiler inline even without @inline?
"Construct a SGKmer from a SGSequence"
function sgkmer(seq::SGSequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(SGKmer{length(seq)}, seq)
end



# Basic Functions
# ---------------

function =={T, K}(a::NucleotideSequence{T}, b::Kmer{T, K})
    if length(a) != K
        return false
    end

    for (u, v) in zip(a, b)
        if u != v
            return false
        end
    end

    return true
end


function =={T, K}(a::Kmer{T, K}, b::NucleotideSequence{T})
    return b == a
end

function Base.getindex{T, K}(x::Kmer{T, K}, i::Integer)
    if i < 1 || i > K
        throw(BoundsError())
    end
    return convert(T, (convert(UInt64, x) >>> (2*(i-1))) & 0b11)
end


function Base.show{K}(io::IO, x::SGKmer{K})
    write(io, "SG $(K)-mer:\n ")
    for i in 1:K
        write(io, convert(Char, x[i]))
    end
end



Base.isless{T, K}(x::Kmer{T, K}, y::Kmer{T, K}) = convert(UInt64, x) < convert(UInt64, y)

Base.length{T, K}(x::Kmer{T, K}) = K

# Iterating over nucleotides
Base.start(x::Kmer) = 1

function Base.next{T, K}(x::Kmer{T, K}, i::Int)
    nt = convert(T, (convert(UInt64, x) >>> (2*(i-1))) & 0b11)
    return (nt, i + 1)
end

Base.done{T, K}(x::Kmer{T, K}, i::Int) = i > K


# Other functions
# ---------------

"""
`complement(kmer::Kmer)`

The Kmer complement of `kmer`
"""
function complement{T, K}(x::Kmer{T, K})
    return convert(Kmer{T, K},
        (~convert(UInt64, x)) & (0xffffffffffffffff >>> (2 * (32 - K))))
end

"""
`reverse(kmer::Kmer)`

Reversed copy of `kmer`
"""
reverse{T, K}(x::Kmer{T, K}) = convert(Kmer{T, K}, nucrev(convert(UInt64, x)) >>> (2 * (32 - K)))

"""
`reverse_complement(kmer::Kmer)`

Reversed complement of `kmer`
"""
reverse_complement{T, K}(x::Kmer{T, K}) = complement(reverse(x))

"""
`mismatches(x::Kmer, y::Kmer)`

Return the number of mismatches between `x` and `y`.

### Arguments
* `a`: first sequence to compare
* `b`: second sequence to compare

### Returns
The number of mismatches
"""
mismatches{T, K}(a::Kmer{T, K}, b::Kmer{T, K}) = nucmismatches(convert(UInt64, a), convert(UInt64, b))

"""
Canonical k-mer of `x`

A canonical k-mer is the numerical lesser of a k-mer and its reverse complement.
This is useful in hashing/counting k-mers in data that is not strand specific,
and thus observing k-mer is equivalent to observing its reverse complement.
"""
function canonical{T, K}(x::Kmer{T, K})
    y = reverse_complement(x)
    return x < y ? x : y
end


"""
Iterate through k-mers neighboring `x` on a de Bruijn graph.
"""
function neighbors{T, K}(x::Kmer{T, K})
    return KmerNeighborIterator{T, K}(x)
end


immutable KmerNeighborIterator{T, K}
    x::Kmer{T, K}
end


Base.start(it::KmerNeighborIterator) = UInt64(0)
Base.done(it::KmerNeighborIterator, i) = i == 4
Base.length(::KmerNeighborIterator) = 4
Base.next{T, K}(it::KmerNeighborIterator{T, K}, i) =
    convert(Kmer{T, K}, (convert(UInt64, it.x) >>> 2) | (i << (2 * K - 2))), i + 1


# EachKmerIterator and EachKmerIteratorState
# ==========================================

# Iterate through every k-mer in a nucleotide sequence
immutable EachKmerIterator{T, K}
    seq::NucleotideSequence{T}
    step::Int
end


immutable EachKmerIteratorState{T, K}
    i::Int
    npos::Int
    x::UInt64
end

# Maybe this function should replace the default constructor.
# Is the (unsafe) default constructor used throughout our code?
"""
Initialize an iterator over all k-mers in a sequence.

Any k-mer containing an N will be skipped over.

### Arguments
  * `t`: Kmer type to enumerate.
  * `seq`: A NucleotideSequence
  * `step`: number of positions between iterated k-mers (default: 1)

### Returns
A EachKmerIterator constructed with these parameters

### Examples
```{.julia execute="false"}
# iterate over codons
for x in each(SGKmer{3}, sg"ATCCTANAGNTACT", 3)
    @show x
end
```
"""
function each{T, K}(t::Type{Kmer{T, K}}, seq::NucleotideSequence{T}, step::Integer=1)
    @assert K >= 0 "K must be ≥ 0 in EachKmer"
    @assert K <= 32 "K must be ≤ 32 in EachKmer"
    @assert step >= 1 "step must be ≥ 1"

    return EachKmerIterator{T, K}(seq, step)
end


function eachkmer{T}(seq::NucleotideSequence{T}, k::Integer, step::Integer=1)
    return each(Kmer{T, K}, seq, step)
end


@inline function nextkmer{T, K}(it::EachKmerIterator{T, K},
                                state::EachKmerIteratorState{T, K}, skip::Int)
    i = state.i + 1
    npos = state.npos
    x = state.x >>> (2 * skip)

    if i > it.seq.part.stop
        return EachKmerIteratorState{T, K}(i, npos, x)
    end

    shift = 2 * (K - 1)
    d, r = divrem32(i - 1)
    data_d = it.seq.data[d + 1] >> 2*r

    lastn = 0
    @inbounds while true
        if i == npos
            npos = nextone(it.seq.ns, i + 1)
            lastn = i
        else
            shift = 2 * (K - skip)
            x |= (((data_d) & 0b11) << shift)
        end

        skip -= 1
        if skip == 0
            if i - lastn < K
                x >>= (2 * it.step)
                skip = it.step
            else
                break
            end
        end

        i += 1
        if i > it.seq.part.stop
            break
        end

        r += 1
        data_d >>= 2
        if r == 32
            r = 0
            d += 1
            if d < length(it.seq.data)
                data_d = it.seq.data[d + 1]
            end
        end
    end

    return EachKmerIteratorState{T, K}(i, npos, x)
end


@inline function start{T, K}(it::EachKmerIterator{T, K})
    i = it.seq.part.start
    npos = nextone(it.seq.ns, i)
    state = EachKmerIteratorState{T, K}(i - 1, npos, UInt64(0))
    if i <= it.seq.part.stop
        return nextkmer(it, state, K)
    else
        return EachKmerIteratorState{T, K}(i, npos, 0)
    end
end


@inline function next{T, K}(it::EachKmerIterator{T, K},
                            state::EachKmerIteratorState{T, K})
    value = convert(Kmer{T, K}, state.x)
    next_state = nextkmer(it, state, it.step)
    return (state.i - it.seq.part.start - K + 2, value), next_state
end


@inline function done{T, K}(it::EachKmerIterator{T, K},
                            state::EachKmerIteratorState{T, K})
    return state.i > it.seq.part.stop
end


# Nucleotide Composition
# ======================

type NucleotideCounts{T <: Nucleotide}
    a::UInt
    c::UInt
    g::UInt
    t::UInt # also hold 'U' count when T == DEPRECATEDNucleotide
    n::UInt

    function NucleotideCounts()
        new(0, 0, 0, 0, 0)
    end
end


# Aliases
# -------

typealias SGNucleotideCounts NucleotideCounts{SGNucleotide}


# Constructors
# ------------

# Count A, C, T/U, G respectively in a kmer stored in a UInt64
function count_a(x::UInt64)
    xinv = ~x
    return count_ones(((xinv >>> 1) & xinv) & 0x5555555555555555)
end
count_c(x::UInt64) = count_ones((((~x) >>> 1) & x) & 0x5555555555555555)
count_g(x::UInt64) = count_ones(((x >>> 1) & (~x)) & 0x5555555555555555)
count_t(x::UInt64) = count_ones((x    & (x >>> 1)) & 0x5555555555555555)

"""
`NucleotideCounts(seq::NucleotideSequence)`

Constructs a NucleotideCounts object from a NucleotideSequence `seq`.
"""
function NucleotideCounts{T}(seq::NucleotideSequence{T})
    dn, rn = divrem64(seq.part.start - 1)

    d = 2*dn
    r = 2*dn

    i = 1
    counts = NucleotideCounts{T}()

    # count leading unaligned bases
    for i in 1:r
        counts[seq[i]] += 1
        i += 1
    end
    if r > 0
        d += 1
    end

    # maybe just skip over blocks of Ns as I go?
    while i + 63 <= length(seq)
        # handle the all-N case
        if seq.ns.chunks[dn + 1] == 0xffffffffffffffff
            counts.n += 64
        else
            counts.a += count_a(seq.data[d + 1]) + count_a(seq.data[d + 2])
            counts.c += count_c(seq.data[d + 1]) + count_c(seq.data[d + 2])
            counts.g += count_g(seq.data[d + 1]) + count_g(seq.data[d + 2])
            counts.t += count_t(seq.data[d + 1]) + count_t(seq.data[d + 2])

            x = seq.ns.chunks[dn + 1]
            if x != 0
                for j in 1:64
                    if x & 0x01 != 0
                        counts.n += 1
                        counts[getnuc(T, seq.data, seq.part.start + i + j - 2)] -= 1
                    end

                    x >>= 1
                    if x == 0
                        break
                    end
                end
            end
        end

        dn += 1
        d += 2
        i += 64
    end

    # count trailing unaligned bases
    while i <= length(seq)
        counts[seq[i]] += 1
        i += 1
    end

    return counts
end

"""
`NucleotideCounts(seq::Kmer)`

Constructs a NucleotideCounts object from a Kmer `seq`.
"""
function NucleotideCounts{T,K}(seq::Kmer{T, K})
    x         = convert(UInt64, seq)
    counts    = NucleotideCounts{T}()
    counts.a += count_a(x) - 32 + K # Take leading zeros into account
    counts.c += count_c(x)
    counts.g += count_g(x)
    counts.t += count_t(x)
    return counts
end



# Basic Functions
# ---------------

Base.setindex!{T}(counts::NucleotideCounts{T}, nt::T) = getfield(counts, Int(convert(UInt, nt) + 1))
Base.setindex!{T}(counts::NucleotideCounts{T}, c::Integer, nt::T) = setfield!(counts, Int(convert(UInt, nt) + 1), c)

# Pad strings so they are right justified when printed
function format_counts(xs)
    strings = AbstractString[string(x) for x in xs]
    len = maximum(map(length, strings))
    for i in 1:length(strings)
        strings[i] = string(repeat(" ", len - length(strings[i])), strings[i])
    end
    return strings
end


# Pretty printing of NucleotideCounts
function Base.show(io::IO, counts::SGNucleotideCounts)
    count_strings = format_counts(
        [counts[SG_A], counts[SG_C], counts[SG_G], counts[SG_T], counts[SG_N]])

    write(io,
        """
        SGNucleotideCounts:
          A => $(count_strings[1])
          C => $(count_strings[2])
          G => $(count_strings[3])
          T => $(count_strings[4])
          N => $(count_strings[5])
        """)
end



# Kmer Composition
# ----------------

"""
Count ocurrences of short (<= 32) k-mers in a sequence.

This method uses a dense table to store counts, requiring O(4^k) memory, so is
not recommended for large k-mer sizes.

### Arguments:
  * 'seq`: A NucleotideSequence
  * `step`: K-mers counted are separated by this many nucleotides (deafult: 1)
"""
immutable KmerCounts{T, K}
    data::Vector{UInt32}

    function KmerCounts(seq::NucleotideSequence{T}, step::Integer=1)
        data = zeros(UInt32, 4^K)
        @inbounds for (_, x) in each(Kmer{T, K}, seq, step)
            data[convert(UInt64, x) + 1] += 1
        end
        return new(data)
    end
end

typealias SGKmerCounts{K} KmerCounts{SGNucleotide, K}

function Base.getindex{T, K}(counts::KmerCounts{T, K}, x::Kmer{T, K})
    @inbounds c = counts.data[convert(UInt64, x) + 1]
    return c
end


function Base.show{T, K}(io::IO, counts::KmerCounts{T, K})
    println(io, (T == SGNucleotide ? "SG" : "DEPRECATED"), "KmerCounts{", K, "}:")
    for x in UInt64(1):UInt64(4^K)
        s = convert(AbstractString, convert(Kmer{T, K}, x - 1))
        println(io, "  ", s, " => ", counts.data[x])
    end
end

