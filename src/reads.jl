
function make_fqparser( filename; forcegzip=false )
   fopen = open( filename, "r" )
   if isgzipped( filename ) || forcegzip
      to_open = ZlibInflateInputStream( fopen, reset_on_end=true )
   else
      to_open = BufferedInputStream( fopen )
   end
   FASTQ.Reader( to_open, fill_ambiguous=DNA_A )
end

# modified for Bio v0.2 with tryread_bool!
@inline function read_chunk!( chunk, parser )
   i = 1
   while i <= length(chunk) && !eof(parser)
      read!( parser, chunk[i] )
      i += 1
   end
   while i <= length(chunk)
      pop!(chunk) # clean up if we are at the end
   end
   parser
end

function allocate_chunk( parser; size=10000 )
  chunk = Vector{eltype(parser)}( size )
  for i in 1:length(chunk)
     chunk[i] = eltype(parser)()
  end
  chunk
end

function allocate_fastq_records( size::Int=10000 )
   chunk = Vector{FASTQRecord}( undef, size )
   for i in 1:length(chunk)
      chunk[i] = FASTQRecord()
   end
   chunk
end

function process_reads!( parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant,
                         multi::MultiMapping{SGAlignSingle}, mod::B;
                         bufsize=150, sam=false, qualoffset=33 ) where B <: BiasModel

   reads  = allocate_fastq_records( bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( stdout )
      write_sam_header( stdbuf, lib )
   end
   while length(reads) > 0
      read_chunk!( reads, parser )
      total += length(reads)
      @inbounds for i in 1:length(reads)
         fill!( reads[i], qualoffset )
         align = ungapped_align( param, lib, reads[i] )
         if !isnull( align )
            biasval = count!( mod, reads[i].sequence )
            if length( align.value ) > 1
               push!( multi, align.value, biasval, quant, lib )
               sam && write_sam( stdbuf, reads[i], align.value, lib, qualoffset=qualoffset )
            else
               count!( quant, align.value[1], biasval )
               sam && write_sam( stdbuf, reads[i], align.value[1], lib, qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(reads[i].sequence) - mean_readlen) / mapped
         end
      end
      if total % 100000 == 0
         GC.gc()
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end


function process_paired_reads!( fwd_parser, rev_parser, param::AlignParam,
                                lib::GraphLib, quant::GraphLibQuant,
                                multi::MultiMapping{SGAlignPaired}, mod::B;
                                bufsize=50, sam=false, qualoffset=33 ) where B <: BiasModel

   fwd_reads  = allocate_fastq_records( bufsize )
   rev_reads  = allocate_fastq_records( bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( stdout )
      write_sam_header( stdbuf, lib )
   end
   while length(fwd_reads) > 0 && length(rev_reads) > 0
      read_chunk!( fwd_reads, fwd_parser )
      read_chunk!( rev_reads, rev_parser )
      total += length(fwd_reads)
      @inbounds for i in 1:length(fwd_reads)
         fill!( fwd_reads[i], qualoffset )
         fill!( rev_reads[i], qualoffset )
         fwd_aln,rev_aln = ungapped_align( param, lib, fwd_reads[i], rev_reads[i] )
         if !isnull( fwd_aln ) && !isnull( rev_aln )
            biasval = count!( mod, fwd_reads[i].sequence, rev_reads[i].sequence )
            if length( fwd_aln.value ) > 1
               push!( multi, fwd_aln.value, rev_aln.value, biasval, quant, lib )
               sam && write_sam( stdbuf, fwd_reads[i], rev_reads[i], fwd_aln.value, rev_aln.value, lib,
                                 paired=true, is_pair_rc=param.is_pair_rc, qualoffset=qualoffset )
            else
               count!( quant, fwd_aln.value[1], rev_aln.value[1], biasval )
               sam && write_sam( stdbuf, fwd_reads[i], fwd_aln.value[1], lib,
                                 paired=true, fwd_mate=true, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
               sam && write_sam( stdbuf, rev_reads[i], rev_aln.value[1], lib,
                                 paired=true, fwd_mate=false, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(fwd_reads[i].sequence) - mean_readlen) / mapped
         end
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end
