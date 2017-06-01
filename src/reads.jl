
function make_fqparser( filename; encoding=Bio.Seq.ILLUMINA18_QUAL_ENCODING, forcegzip=false )
   fopen = open( filename, "r" )
   if isgzipped( filename ) || forcegzip
      to_open = ZlibInflateInputStream( fopen, reset_on_end=true )
   else
      to_open = BufferedInputStream( fopen )
   end 
   FASTQReader{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{2}}}( to_open, 
                                                             encoding,
                                                             DNA_A ), Requests.ResponseStream{TCPSocket}()
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

# This function specifically tries to download a fastq file from a url string and returns
# the Bio.Seq parser, the IOBuffer, and the RemoteChannel to the HTTPC.get
function make_http_fqparser( url::String; encoding=Bio.Seq.ILLUMINA18_QUAL_ENCODING, forcegzip=false )
   response = Requests.get_streaming(url)
   if isgzipped( url ) || forcegzip
      zlibstr  = ZlibInflateInputStream( response.buffer, reset_on_end=true )
      fqparser = FASTQReader{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{2}}}( zlibstr,         encoding, DNA_A )
   else
      fqparser = FASTQReader{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{2}}}( response.buffer, encoding, DNA_A )
   end
   fqparser, response
end

# Use this version to parse reads from a parser that is reliant on the state
function read_http_chunk!( chunk, parser, resp; maxtime=24 )
   i = 1
   const iobuf      = resp.buffer
   const nb_needed  = 8192
   const start_mark = iobuf.mark
   const start_size = iobuf.size
   const start_time = time()
   if !(200 <= resp.response.status < 300)
      error("HTTP Code $(resp.response.status)! Download failed!")
   end
   while i <= length(chunk) && !(eof(iobuf) && eof(parser))
      if resp.state!=Requests.BodyDone && nb_available(iobuf) < nb_needed
         sleep(eps(Float64))
         if time() - start_time > maxtime
            error("HTTP Timeout! Unable to download file!")
         end
         continue
      end
      read!( parser, chunk[i] )
      i += 1
   end
   while i <= length(chunk)
      pop!(chunk) # clean up if we are at the end
   end
   parser
end

function allocate_chunk( parser; size=100000 )
  chunk = Vector{eltype(parser)}( size )
  for i in 1:length(chunk)
     chunk[i] = eltype(parser)()
  end
  chunk
end


function process_reads!( parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant, 
                         multi::Vector{Multimap}; bufsize=50, sam=false, qualoffset=33,
                         response=Requests.ResponseStream{TCPSocket}(), http=false )
  
   const reads  = allocate_chunk( parser, size=bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( STDOUT )
      write_sam_header( stdbuf, lib )
   end
   while length(reads) > 0
      if http
         read_http_chunk!( reads, parser, response )
      else
         read_chunk!( reads, parser )
      end
      total += length(reads)
      for i in 1:length(reads)
         align = ungapped_align( param, lib, reads[i] )
         if !isnull( align )
            if length( align.value ) > 1
               push!( multi, Multimap( align.value ) )
               sam && write_sam( stdbuf, reads[i], align.value, lib, qualoffset=qualoffset )
            else
               count!( quant, align.value[1] )
               sam && write_sam( stdbuf, reads[i], align.value[1], lib, qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(reads[i].seq) - mean_readlen) / mapped
         end
      end
      if total % 100000 == 0
         gc()
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end


function process_paired_reads!( fwd_parser, rev_parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant,
                                multi::Vector{Multimap}; bufsize=50, sam=false, qualoffset=33,
                                     response=Requests.ResponseStream{TCPSocket}(), 
                                mate_response=Requests.ResponseStream{TCPSocket}(), 
                                http=false )

   const fwd_reads  = allocate_chunk( fwd_parser, size=bufsize )
   const rev_reads  = allocate_chunk( rev_parser, size=bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( STDOUT )
      write_sam_header( stdbuf, lib )
   end
   while length(fwd_reads) > 0 && length(rev_reads) > 0
      if http
         read_http_chunk!( fwd_reads, fwd_parser, response )
         read_http_chunk!( rev_reads, rev_parser, mate_response )
      else
         read_chunk!( fwd_reads, fwd_parser )
         read_chunk!( rev_reads, rev_parser )
      end
      total += length(fwd_reads)
      for i in 1:length(fwd_reads)
         fwd_aln,rev_aln = ungapped_align( param, lib, fwd_reads[i], rev_reads[i] )
         if !isnull( fwd_aln ) && !isnull( rev_aln )
            if length( fwd_aln.value ) > 1
               push!( multi, Multimap( fwd_aln.value ) )
               push!( multi, Multimap( rev_aln.value ) )
               sam && write_sam( stdbuf, fwd_reads[i], rev_reads[i], fwd_aln.value, rev_aln.value, lib,
                                 paired=true, is_pair_rc=param.is_pair_rc, qualoffset=qualoffset )
            else
               count!( quant, fwd_aln.value[1], rev_aln.value[1] )
               sam && write_sam( stdbuf, fwd_reads[i], fwd_aln.value[1], lib, 
                                 paired=true, fwd_mate=true, is_pair_rc=param.is_pair_rc, 
                                 qualoffset=qualoffset )
               sam && write_sam( stdbuf, rev_reads[i], rev_aln.value[1], lib, 
                                 paired=true, fwd_mate=false, is_pair_rc=param.is_pair_rc, 
                                 qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(fwd_reads[i].seq) - mean_readlen) / mapped
         end
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end
