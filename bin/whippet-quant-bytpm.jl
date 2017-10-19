#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = readline(open(dir * "/VERSION"))

tic()
println( STDERR, "Whippet $ver loading and compiling... " )

using ArgParse

unshift!( LOAD_PATH, dir * "/../src" )
using Whippet

function parse_cmd()
   s = ArgParseSettings()
   @add_arg_table s begin
    "--tpm"
      help     = "File containing two columns, isoform names, and transcripts-per-million (TpM)"
      arg_type = String
    "--index", "-x"
      help     = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph[.jls])"
      arg_type = String
      default  = fixpath( "$dir/../index/graph" )
    "--out", "-o"
      help     = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = String
      default  = fixpath( "./output" )
   end
   return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, " $( round( toq(), 6 ) ) seconds" )

   indexpath = fixpath( args["index"] )
   println(STDERR, "Loading splice graph index... $( indexpath ).jls")
   @timer const lib = open(deserialize, "$( indexpath ).jls")

   tpms = Vector{Vector{Float64}}(length(lib.graphs))
   fill_tpms!( tpms )

   println(STDERR, "Whippet $ver done." )
end

function fill_tpms!( , )
   
end

@timer main()
