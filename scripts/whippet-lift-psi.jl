#!/usr/bin/env julia
# Tim Sterne-Weiler 2016

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/../bin/VERSION")))

tic()
println( STDERR, "Whippet $ver loading and compiling... " )

using Bio.Intervals
using IntervalTrees
using Libz
using BufferedStreams
using ArgParse
using Glob

include("../src/types.jl")
include("src/bio_ext.jl")
include("src/chain.jl")

function parse_cmd()
  s = ArgParseSettings()
  # TODO finish options...
  @add_arg_table s begin
    "--to"
      help     = ".psi.gz file to liftover to"
      arg_type = String
      required = true
    "--from"
      help     = ".psi.gz file to liftover"
      arg_type = String
      required = true
    "--chain"
      help     = ".chain.gz file to use for liftover"
      arg_type = String
      required = true
    "--out", "-o"
      help     = "Core file name to send .diff.gz output to!"
      arg_type = String
      default  = fixpath( "./output" )
  end
  return parse_args(s)
end

function main()
   args  = parse_cmd()
   println(STDERR, "Converting coordinates...") 
   
   chain = load_io_constructor( args["chain"], LiftOverChain )
   from  = load_io_constructor( args["from"],  PsiCollection )
   to    = load_io_constructor( args["to"],    PsiCollection )

   println(STDERR, "Whippet $ver done." )   
end

function load_io_constructor( filename::String, constructor::Type )
    os = BufferedInputStream( open( fixpath(filename) ) )
    if isgzipped( filename )
        stream = ZlibInflateInputStream( os )
    else
        stream = os
    end
    chain = constructor( stream )
    close(stream)
    close(os)
    chain
end

immutable PsiEntry{S <: AbstractString}
   gene::S
   node::Int64
   event::S
   Psi::Float64
   ci_width::Float64
   ci_range::S # TODO make tuple->float,float
   complexity::S
   entropy::Float64
   inc::S
   exc::S
end

typealias PsiCollection Bio.Intervals.IntervalCollection{PsiEntry}

function PsiCollection( io )

   psi = PsiCollection()

   for l in eachline( io )
      spl      = split(chomp(l), '\t')
      spl[7] == "NA" && continue
      entry    = PsiEntry( spl[1], parse(Int, spl[2]), spl[5], 
                           parse(Float64, spl[6]),
                           parse(Float64, spl[7]),
                           spl[8],spl[9],
                           parse(Float64, spl[10]),
                           spl[11], spl[12] )
      psi_line = Bio.Intervals.Interval(String(spl[3]), spl[4][1], entry)
      push!( psi, psi_line )
   end
   psi
end

main()
