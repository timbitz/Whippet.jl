#!/usr/bin/env julia
# Tim Sterne-Weiler 2018

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "../VERSION")))

tic()
println( STDERR, "Whippet $ver loading and compiling... " )

using ArgParse

push!( LOAD_PATH, dir * "/../../src" )
using Whippet

using BufferedStreams
using Libz

include("chain.jl")

function parse_cmd()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--psi"
      help     = "Psi input file"
      arg_type = String
    "--chain"
      help     = "Chain file"
      arg_type = String
  end
  return parse_args(s)
end

function main()
   args  = parse_cmd()

   psiio = open( args["psi"] )
   psibuf = isgzipped( args["psi"] ) ? ZlibInflateInputStream( psiio ) : BufferedInputStream( psiio )

   chainio = open( args["chain"] )
   chainbuf = isgzipped( args["psi"] ) ? ZlibInflateInputStream( chainio ) : BufferedInputStream( chainio )

   lift = LiftOverChain( chainbuf )

   # iterate through psifile
   for ln in eachline( psibuf )
      # collect each node into Interval( chr, start, end, strand )
      # push! to interval array
   end
   # liftover interval array

   println(STDERR, "Whippet $ver done." )
end

@timer main()
