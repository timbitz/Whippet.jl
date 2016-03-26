
tab_write{S <: AbstractString}( io::BufOut, str::S ) = (write( io, str ); write( io, '\t'  ))
tab_write( io::BufOut, str::Char ) = (write( io, str ); write( io, '\t'  ))

function coord_write( io::BufOut, chr, first, last )
   write( io, chr   )
   write( io, ':'   )
   write( io, string(first) )
   write( io, '-'   )
   write( io, string(last)  )
end

function coord_write( io::BufOut, chr, first, last, strand )
   coord_write( io, chr, first, last )
   write( io, ':'    )
   write( io, strand )
end

function sam_flag( align::SGAlignment, lib::GraphLib, ind )
   flag = UInt16(0)
   lib.info[ ind ].strand || flag &= 0x10
   
end

function cigar_string( align::SGAlignment, sg::SpliceGraph )
   matches = align.matches
   cigar = ""
   for i in 1:length( align.path )-1
      mat = min( matches + align.offset, sg.offset[i] + sg.nodelen[i] )
      matches -= mat
      cigar *= string(mat) * 'M'
      const intron = sg.offset[i+1] - sg.offset[i] + sg.nodelen[i] - 2
      if intron > 0
         cigar *= string(intron) * 'N'
      end
   end   
   
end

function write_sam( io::BufOut, read::SeqRecord, align::SGAlignment, lib::GraphLib, )
   const geneind = align.path[1].gene
   const sg = lib.graphs[geneind] 
   tab_write( io, read.name )
   tab_write( io, string( sam_flag(align, lib, geneind) ) )
   tab_write( io, lib.info[geneind].name )
   tab_write( io, string(lib.offset[geneind] + sg.offset[align.path[1].node] + align.offset) )
   tab_write( io, '0' )
   tab_write( io, cigar_string( align, sg ) )
end

