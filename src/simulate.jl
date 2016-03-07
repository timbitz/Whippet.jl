#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using ArgParse

dir = splitdir(@__FILE__)[1]

push!( LOAD_PATH, dir * "/../src" )
import SpliceGraphs
@everywhere using SpliceGraphs

function parse_cmd()
  s = ArgParseSettings(version="0.0.1", add_version=true)
  # TODO finish options...
  @add_arg_table s begin
    "--index", "-x"
      help = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = ASCIIString
      default  = fixpath( "$(dir)/../index/graph" )
    "--out", "-o"
      help = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = ASCIIString
      default  = fixpath( "$(dir)/../output" )
    "--max-complexity", "-c"
      help = "What is the maximum complexity we should allow?"
      arg_type = Int64
      default  = 8
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, "Loading splice graph index... $( args["index"] ).jls")
   @time const lib = open(deserialize, "$( args["index"] ).jls")

   println(STDERR, "Loading annotation index... $( args["index"] )_anno.jls")
   @time const anno = open(deserialize, "$( args["index"] )_anno.jls")

end

function simulate_genes( lib, anno, max_comp )
   for g in 1:length(lib.graphs)
      simulate_psi( lib.graphs[g], lib.names[g], max_comp )
   end
end

function simulate_psi( sg::SpliceGraph, gene::ASCIIString, comp::Int )
   max_comp = min( comp, length(sg.nodelen) - 2, 0 )
   for n in 1:length(sg.nodelen)
     

   end
end

main()
