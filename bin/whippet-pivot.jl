#!/usr/bin/env julia
# Tim Sterne-Weiler 2021

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

start = time_ns()
println( stderr, "Whippet $ver loading... " )

Pkg.activate(dir * "/..")

using Whippet
using Random
using ArgParse
using Glob


function parse_cmd()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--files", "-f"
      help     = "Pattern to glob.psi[.gz] (common-filename-segment [*.psi*]), or comma delimited list of filenames. ie. (-a sample_a) would work for sample_a-rep1.psi.gz,sample_a-rep2.psi.gz,..."
      arg_type = String
      default  = ""
#      required = true
    "--directory", "-d"
      help     = "Directory to search for file patterns or list in -a and -b"
      arg_type = String
      default  = "."
  end
  return parse_args(s)
end

function retrievefilelist( pattern::String, dir::String )
   list = Vector{String}()
   if something(findfirst(isequal(','), pattern), 0) > 0
      tmp = split( pattern, ',', keepempty=false )
   else
        pat = length(pattern) > 0 ? "*" * pattern : pattern
      tmp = glob( pat * "*.psi*", dir )
   end
   # now clean the return
   for file in tmp
      push!( list, string(file) )
   end
   list
end

readlinesplit( stream ) = readline(stream) |> x->split( x, '\t' )

function gene_centric( streams::Vector{BufferedStreams.BufferedInputStream},
                       names::Vector{String} )
   
   curline  = Vector{Vector{SubString{String}}}(undef, length(streams))
   metadata = Vector{String}()
   header   = ""

   # initialize header lines
   for (i,s) in enumerate(streams)
   	header = readline(s)
   	curline[i] = readlinesplit(s) # first psi line
   end
   curgene = curline[1][1]
   ostream = ZlibDeflateOutputStream(open(curgene * ".gpsi", "w"))
   write(ostream, "Sample\t" * header)

   while !eof(streams[1])
   	for (i,s) in enumerate(streams)
      	while curline[i][1] == curgene && !eof(s)
      		tab_write(ostream, names[i])
      		tab_write(ostream, curline[i])

      		curline[i] = readlinesplit(s)
      	end
      end
      close(ostream)
      curgene = curline[1][1]
   	ostream = ZlibDeflateOutputStream(open(curgene * ".gpsi", "w"))
      write(ostream, "Sample\t" * header)
      println(stderr, "Pivoting into " * curgene * ".gpsi")
   end
end

function main()
   args  = parse_cmd()
   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )

   dir   = fixpath( args["directory"] )
   files = retrievefilelist( args["files"], dir )
   println(stderr, "Loading files to pivot: $(join(map(basename, files), '\n'))")

   if length(files) <= 0
      error("Unable to find any files! length(files) == $(length(files))")
   end

   fstreams = open_streams( files )

   println(stderr, "Pivoting to gene-centric psi files...")
 
   gene_centric(fstreams)

   println(stderr, "Whippet $ver done." )
end

@timer main()