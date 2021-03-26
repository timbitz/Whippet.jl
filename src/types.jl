#=
Basic Types and Aliases
=#

const Mb = 1_000_000
const GENOMESIZE = 3235Mb

# ALIAS NEW TYPES FOR INCREASED CODE READABILITY
if GENOMESIZE < typemax(UInt32)
   const CoordInt = UInt32
else
   const CoordInt = UInt64
end

const CoordTuple = Tuple{Vararg{CoordInt}}
const CoordArray = Vector{CoordInt}
const CoordTree  = IntervalTree{CoordInt,IntervalTrees.Interval{CoordInt}}
const CoordSet   = Set{Tuple{CoordInt,CoordInt}}

const ExonInt    = UInt16
const ExonTree   = IntervalTree{ExonInt,IntervalTrees.Interval{ExonInt}}

const SGSequence = BioSequences.BioSequence{DNAAlphabet{4}}

const GeneName    = String
const SeqName     = String
const RefSeqId    = String

const GeneMeta = Tuple{GeneName, SeqName, Char}

const BufOut = BufferedStreams.BufferedOutputStream

struct GeneInfo
   gene::GeneName
   name::SeqName
   strand::Bool # pos is true, neg is false
end

GeneInfo( gene::S, name::S, strand::Char ) where {S <: AbstractString} =
                               GeneInfo( convert(GeneName, gene),
                                         convert(SeqName, name),
                                         strand == '+' ? true :
                                                         false )

struct TxInfo
   name::RefSeqId
   txstart::CoordInt
   txend::CoordInt
   exnum::CoordInt
end


#=
Basic IO Functions
=#

fixpath( str::String ) = abspath( expanduser( str ) )

isgzipped( filename::String ) = splitext(filename)[2] == ".gz"
hasextension( filename::String, ext::String ) = splitext(filename)[2] == ext

function increment!( dict::Dict{K,V}, key::K, val::V=one(V) ) where {K,V}
   if haskey( dict, key )
      dict[key] += val
   else
      dict[key] = val
   end
   dict
end

# add constructor for compatibility
BioSequence{A}()             where {A <: Alphabet} = LongSequence{A}()
BioSequence{A}( var )        where {A <: Alphabet} = LongSequence{A}(var)
BioSequence{A}( arr, a, b )  where {A <: Alphabet} = LongSequence{A}(arr, a, b)
