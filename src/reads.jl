
# TODO fix type instability
function reduce{T <: Vector{Nullable}}( a::T, b::T )
   if isnull(a)
      return b
   elseif isnull(b)
      return a
   else
      return [a;b]
   end
end

function make_fqparser( dir, filename )
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

function allocate_chunk( parser, size=200_000 )
  chunk = Vector{eltype(parser)}( size )
  for i in 1:length(chunk)
     chunk[i] = eltype(parser)()
  end
  chunk
end

function read_chunk!( chunk, parser )
   i = 1
   while i <= length(chunk) && read!( parser, chunk[i] )
      i += 1
   end
   while i < length(chunk)
      pop!(chunk) # clean up if we are at the end
   end
end
