function chunk_ranges( datasize, num=nworkers() )
   size,left = divrem( datasize, num )
   len = fill( size, num )
   len[1:left] += 1
   offset = 1
   ranges = UnitRange[]
   for i in 1:num
      push!(ranges, offset:(offset+len[i]-1))
      offset += len[i]
   end
   ranges
end

function process_reads!( parser, param::AlignParam, lib::GraphLib,
                         quant::GraphLibQuant, multi::Vector{Multimap}; bufsize=100)
   if VERSION >= v"0.5.0-dev" && nthreads() > 1
      return _process_reads_faster!( parser, param, lib, quant, multi, bufsize=bufsize)
   else
      return _process_reads!( parser, param, lib, quant, multi, bufsize=bufsize)
   end
end

function _process_reads_faster!( parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant,
                                 multi::Vector{Multimap}; bufsize=100 )
   reads  = allocate_chunk( parser, bufsize )
   align  = Vector{SGAlignVec}( length(reads) )
   #ranges = chunk_ranges( bufsize, nthreads() )
   mean_readlen = 0.0
   total  = 0
   mapped = 0
   while length(reads) > 0
      read_chunk!( reads, parser )
      @threads for i in 1:length(reads)
         align[i] = ungapped_align( param, lib, reads[i] )
      end
      for i in 1:length(reads)
         if !isnull( align[i] )
            if length( get(align[i]) ) > 1
               push!( multi, Multimap( get(align[i]) ) )
            else
               count!( quant, get(align[i])[1] )
            end
            mapped += 1
         end
         total += 1
         @fastmath mean_readlen += (length(reads[i].seq) - mean_readlen) / total
      end
   end # end while
   mapped,total,mean_readlen
end

