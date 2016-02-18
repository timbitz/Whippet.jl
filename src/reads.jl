
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

function allocate_chunk( parser, size=250000 )
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

function process_reads!( parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant, 
                         multi::Vector{Multimap}; bufsize=100000 )
   mapped = 0
   unmapped = 0
   reads  = allocate_chunk( parser, bufsize )
   while length(reads) > 0
      read_chunk!( reads, parser )
      for r in reads
         align = ungapped_align( param, lib, r )
         if !isnull( align )
            if length( get(align) ) > 1
               push!( multi, Multimap( get(align) ) )
            else
               count!( quant, get(align) )
            end
            mapped += 1
         else
            unmapped += 1
         end
      end
      nreads += length(reads)
   end # end while
   mapped,unmapped
end
