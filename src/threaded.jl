
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
   mapped = 0
   reads  = allocate_chunk( parser, bufsize )
   mean_readlen = 0.0
   total = 0
   mut = Mutex()
   while length(reads) > 0
      read_chunk!( reads, parser )
      @threads for i in 1:length(reads)
         local align = ungapped_align( param, lib, reads[i] )
         lock!(mut)
         if !isnull( align )
            if length( get(align) ) > 1
               push!( multi, Multimap( get(align) ) )
            else
               count!( quant, get(align)[1] )
            end
            mapped += 1
         end
         total += 1
         unlock!(mut)
         @fastmath mean_readlen += (length(reads[i].seq) - mean_readlen) / total
      end
   end # end while
   mapped,total,mean_readlen
end

