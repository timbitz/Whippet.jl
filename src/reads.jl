
function make_fqparser( filename )
   re = match( r"(\S+).(fq|fastq)(.gz?)", filename )
   if re != nothing
      if re.captures[3] != nothing #then gzipped
         to_open = open( filename ) |> ZlibInflateInputStream
      else
         to_open = filename
      end
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


