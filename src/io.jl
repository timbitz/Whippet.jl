
tab_write{S <: AbstractString}( io::BufOut, str::S ) = (write( io, str ); write( io, '\t'  ))
tab_write( io::BufOut, str::Char ) = (write( io, str ); write( io, '\t'  ))

function coord_write( io::BufOut, chr, first, last; tab=false )
   write( io, chr   )
   write( io, ':'   )
   write( io, string(first) )
   write( io, '-'   )
   write( io, string(last)  )
   tab && write( io, '\t' )
end

function coord_write( io::BufOut, chr, first, last, strand; tab=false )
   coord_write( io, chr, first, last )
   write( io, ':'    )
   write( io, strand )
   tab && write( io, '\t' )
end

function seq_write( io::BufOut, read::SeqRecord; tab=false )
   for c in read.seq
      write( io, c )
   end
   tab && write( io, '\t' )
end

function qual_write( io::BufOut, read::SeqRecord, qualoffset=64; tab=false )
   for i in read.metadata.quality
      write( io, convert(Char, i+qualoffset) )
   end
   tab && write( io, '\t' )
end

function sam_flag( align::SGAlignment, lib::GraphLib, ind )
   flag = UInt16(0)
   lib.info[ ind ].strand || (flag &= 0x10)
   
end

function cigar_string( align::SGAlignment, sg::SpliceGraph )
   matchleft = align.matches
   cigar = ""
   curpos = align.offset
   for i in 1:length( align.path )
      if matchleft + curpos <= sg.offset[i] + sg.nodelen[i] - 1
         cigar *= string( matchleft ) * 'M'
      else
         curspace = sg.offset[i] + sg.nodelen[i] - 1 - curpos
         cigar *= string( curspace ) * 'M'
         matchleft -= curspace
         if i < length( align.path )
            const intron = sg.offset[i+1] - sg.offset[i] + sg.nodelen[i] - 1
            if intron > 0
               cigar *= string( intron ) * 'N'
            end
         else
            error("wouldn't be a valid alignment, can't make SAM entry!")
         end
      end
   end
   cigar
end

function write_sam( io::BufOut, read::SeqRecord, align::SGAlignment, lib::GraphLib; mapq=0 )
   const geneind = align.path[1].gene
   const sg = lib.graphs[geneind] 
   tab_write( io, read.name )
   tab_write( io, string( sam_flag(align, lib, geneind) ) )
   tab_write( io, lib.info[geneind].name )
   tab_write( io, string(lib.offset[geneind] + sg.offset[align.path[1].node] + align.offset - 1) )
   tab_write( io, string(mapq) )
   tab_write( io, cigar_string( align, sg ) )
   tab_write( io, '*' )
   tab_write( io, '0' )
   tab_write( io, '0' )
   seq_write( io, read, tab=true )
   qual_write( io, read )
end

