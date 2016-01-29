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

typealias Genename    ASCIIString

