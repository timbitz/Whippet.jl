const Mb = 1_000_000
const GENOMESIZE = 3235Mb

# ALIAS NEW TYPES FOR INCREASED CODE READABILITY
if GENOMESIZE < typemax(UInt32)
   typealias CoordInt UInt32
else
   typealias CoordInt UInt64
end

typealias CoordTuple Tuple{Vararg{CoordInt}}
typealias CoordArray Vector{CoordInt}
typealias CoordTree IntervalTree{CoordInt,Interval{CoordInt}}

typealias ExonMax    UInt16

if VERSION < v"0.5.0-dev"
   typealias String ASCIIString
end
typealias GeneName   String
typealias SeqName    String
typealias RefSeqId   String

typealias GeneMeta   Tuple{GeneName, SeqName, Char}

typealias BufOut BufferedStreams.BufferedOutputStream

immutable GeneInfo
   name::SeqName
   strand::Bool # pos is true, neg is false
end

GeneInfo( name::SeqName, strand::Char ) = GeneInfo( name, strand == '+' ? true : false )

immutable TxInfo
   name::RefSeqId
   txstart::CoordInt
   txend::CoordInt
   exnum::CoordInt
end

