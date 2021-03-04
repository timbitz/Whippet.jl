#!/usr/bin/env julia
# Tim Sterne-Weiler 2016

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
  # TODO finish options...
  @add_arg_table s begin
    "--a", "-a"
      help     = "Replicates for Set A -- Could be: pattern to glob.psi (common-filename-segment [*.psi*]), or comma delimited list of filenames. ie. (-a sample_a) would work for sample_a-rep1.psi.gz,sample_a-rep2.psi.gz,..."
      arg_type = String
#      required = true
    "--b", "-b"
      help     = "Replicates for Set B -- Same rules as for (-a)"
      arg_type = String
#      required = true
    "--out", "-o"
      help     = "Core file name to send .diff.gz output to!"
      arg_type = String
      default  = fixpath( "./output" )
    "--directory", "-d"
      help     = "Directory to search for file patterns or list in -a and -b"
      arg_type = String
      default  = "."
    "--min-reads", "-r"
      help     = "Minimum number of reads for a single event to be included!"
      arg_type = Int64
      default  = 5
    "--min-samples", "-s"
      help     = "Minimum number of samples in a or in b for each event to be considered!"
      arg_type = Int64
      default  = 1
#=    "--min-delta-psi", "-m"
      help     = "Calculate max probability of |deltaPsi| greater than this value (default is 0.0, it is not advisable to change this)."
      arg_type = Float64
      default  = 0.0 =#
    "--emperical-size", "-e"
      help     = "Emperical distribution size to sample from."
      arg_type = Int64
      default  = 1000
    "--seed", "-g"
      help     = "Seed the RNG (Int) for reproducible results on successive runs"
      arg_type = Int64
      default  = 123456
  end
  return parse_args(s)
end

function retrievefilelist( pattern::String, dir::String )
   list = Vector{String}()
   if something(findfirst(isequal(','), pattern), 0) > 0
      tmp = split( pattern, ',', keepempty=false )
   else
      tmp = glob( "*" * pattern * "*.psi*", dir )
   end
   # now clean the return
   for file in tmp
      push!( list, string(file) )
   end
   list
end

function main()
   args  = parse_cmd()
   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )
   Random.seed!( args["seed"] )
   dir   = fixpath( args["directory"] )
   lista = retrievefilelist( args["a"], dir )
   listb = retrievefilelist( args["b"], dir )
   println(stderr, "Sample A: $(join(map(basename, lista), ','))")
   println(stderr, "Sample B: $(join(map(basename, listb), ','))")
   if length(lista) <= 0 || length(listb) <= 0
      error("Unable to match files! length(a) == $(length(lista)), length(b) == $(length(listb))!")
   end

   astreams = open_streams( lista )
   bstreams = open_streams( listb )

   println(stderr, "Now processing files and calculating posterior distributions...")
   @timer process_psi_files( args["out"] * ".diff.gz", astreams, bstreams,
                             min_samp=args["min-samples"],
                             min_reads=args["min-reads"],
                             amt=0.0,
                             size=args["emperical-size"] )
   println(stderr, "Whippet $ver done." )
end

@timer main()
