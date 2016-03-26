const Mb = 1_000_000
const GENOMESIZE = 3235Mb

# ALIAS NEW TYPES FOR INCREASED CODE READABILITY
if GENOMESIZE < typemax(UInt32)
   typealias Coordint UInt32
else
   typealias Coordint UInt64
end
typealias Coordtuple Tuple{Vararg{Coordint}}
typealias Coordarray Vector{Coordint}
typealias Exonmax    UInt16
typealias Genename   ASCIIString
typealias Seqname    ASCIIString
typealias GeneMeta Tuple{Genename, Seqname, Char}

typealias BufOut BufferedStreams.BufferedOutputStream

immutable GeneInfo
   name::Seqname
   strand::Bool # pos is true, neg is false
end

GeneInfo( name::Seqname, strand::Char ) = GeneInfo( name, strand == '+' ? true : false )
