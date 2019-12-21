
mutable struct FASTQRecord
   sequence::BioSequence{DNAAlphabet{2}}
   quality::Vector{UInt8}
   raw::FASTQ.Record

   FASTQRecord() = new( BioSequence{DNAAlphabet{2}}(),
                        Vector{UInt8}(), FASTQ.Record() )
end

Base.read!( parser::FASTQ.Reader, rec::FASTQRecord ) = read!( parser, rec.raw )

function Base.fill!( rec::FASTQRecord, offset=33 )
   rec.sequence = BioSequences.BioSequence{BioSequences.DNAAlphabet{2}}( rec.raw.data,
                                                                         first(rec.raw.sequence),
                                                                         last(rec.raw.sequence) )
   rec.quality  = rec.raw.data[rec.raw.quality] .- convert(UInt8, offset)
   rec
end

Base.string( rec::FASTQRecord ) = String(rec.raw.data[rec.raw.identifier])

function reverse_complement!(rec::FASTQRecord)
   BioSequences.reverse_complement!(rec.sequence)
   reverse!(rec.quality)
end
