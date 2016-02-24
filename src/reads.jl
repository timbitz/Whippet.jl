
function make_fqparser( filename )
   if isgzipped( filename )
      to_open = open( filename, "r" ) |> ZlibInflateInputStream
   else
      to_open = filename
   end 
   open( to_open, FASTQ )
end

function read_chunk!( chunk, parser )
   i = 1
   while i <= length(chunk) && read!( parser, chunk[i] )
      i += 1
   end
   while i <= length(chunk)
      pop!(chunk) # clean up if we are at the end
   end
   parser
end

function allocate_chunk( parser, size=100000 )
  chunk = Vector{eltype(parser)}( size )
  for i in 1:length(chunk)
     chunk[i] = eltype(parser)()
  end
  chunk
end

allocate_rref( size=250000; rreftype=RemoteRef{Channel{Any}} ) = Vector{rreftype}( size )

function resize_rref!( rref, subsize )
   @assert subsize <= length(rref)
   while subsize <= length(rref)
      pop!(rref) # clean up if we are at the end
   end
end

function sendto(p::Int, nm, val)
   ref = @spawnat(p, eval(Main, Expr(:(=), nm, val)))
end

macro sendto(p, nm, val)
   return :( sendto($p, $nm, $val) )
end

macro broadcast(nm, val)
   quote
      @sync for p in workers()
         @async sendto(p, $nm, $val)
      end
   end
end


process_reads!( parser, param::AlignParam, lib::GraphLib,
                quant::GraphLibQuant, multi::Vector{Multimap}; bufsize=100) = _process_reads!( parser, param, lib, quant,
                                                                                               multi, bufsize=bufsize )

function _process_reads!( parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant, 
                         multi::Vector{Multimap}; bufsize=100 )
   mapped = 0
   reads  = allocate_chunk( parser, bufsize )
   mean_readlen = 0.0
   total = 0
   while length(reads) > 0
      read_chunk!( reads, parser )
      for i in 1:length(reads)
         align = ungapped_align( param, lib, reads[i] )
         if !isnull( align )
            if length( get(align) ) > 1
               push!( multi, Multimap( get(align) ) )
            else
               count!( quant, get(align)[1] )
            end
            mapped += 1
         end
      
         total += 1
         @fastmath mean_readlen += (length(reads[i].seq) - mean_readlen) / total
      end
   end # end while
   mapped,total,mean_readlen
end

