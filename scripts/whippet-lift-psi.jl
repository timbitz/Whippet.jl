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
    "--from-meta"
      help     = "from species"
      arg_type = String
      default  = "hg19"
    "--to-meta"
      help     = "to species"
      arg_type = String
      default  = "query"
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
   
   chain  = load_io_constructor( args["chain"], LiftOverChain )
   from   = load_io_constructor( args["from"],  PsiCollection )
   to     = load_io_constructor( args["to"],    PsiCollection )

   from_meta = split_comma( args["from-meta"] )
   to_meta   = split_comma( args["to-meta"] )

   lifted  = liftover( chain, from )
   nonnull = lifted_intervals( lifted )
   println(STDERR, "Intersecting $(length(nonnull)) lifted coordinates from $(length(lifted)) original" )
   for (f,t) in intersect( nonnull, to )
      ref_id = "$(f.metadata.gene):$(f.metadata.node)"
      tmp_from_meta = copy(from_meta)
      tmp_to_meta   = copy(to_meta)
      push!(tmp_from_meta, ref_id)
      push!(tmp_to_meta,   ref_id)
      fcoord = "$(f.seqname):$(f.first)-$(f.last)"
      tcoord = "$(t.seqname):$(t.first)-$(t.last)"
      writepsi(STDOUT, f, tmp_from_meta, fcoord )
      writepsi(STDOUT, t, tmp_to_meta, tcoord )
   end

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


type PsiEntry
   gene::String
   node::Int64
   event::String
   psi::Float64
   ci_width::Float64
   ci_range::String # TODO make tuple->float,float
   total::Float64
   complexity::String
   entropy::Float64
   inc::String
   exc::String
end

typealias PsiCollection Bio.Intervals.IntervalCollection{PsiEntry}

function PsiCollection( io )

   psi = PsiCollection()
   header = true

   for l in eachline( io )
      spl      = split(chomp(l), '\t')
      isna = spl[6] == "NA" ? true : false
      (spl[5] == "BS" || spl[7] == "NA") && continue
      if header
         header = false
         continue
      end
      try
         if !isna
            entry = PsiEntry( String(spl[1]), parse(Int, spl[2]), String(spl[5]), 
                              parse(Float64, spl[6]),
                              parse(Float64, spl[7]),
                              String(spl[8]), parse(Float64, spl[9]), String(spl[10]),
                              parse(Float64, spl[11]),
                              String(spl[12]), String(spl[13]) )
         else
            entry = PsiEntry( String(spl[1]), parse(Int, spl[2]), String(spl[5]),
                              parse(Float64, spl[6]),
                              0.0,
                              String(spl[8]), 0.0, String(spl[10]),
                              0.0,
                              String(spl[12]), String(spl[13]) )
         end
         psi_line = Bio.Intervals.Interval(String(spl[3]), spl[4][1], entry)
         push!( psi, psi_line )
      catch e
         println(STDERR, l)
         rethrow(e)
      end
   end
   psi
end

split_comma{S <: AbstractString}( str::S ) = map( String, split(str, ',', keep=false) )

plug_write{S <: AbstractString}( io, str::S; plug::Char='\t' ) = (write( io, str ); write( io, plug  ))
plug_write( io, str::Char; plug::Char='\t' ) = (write( io, str ); write( io, plug  ))

tab_write{S <: AbstractString}( io, str::S ) = plug_write( io, str, plug='\t' )
tab_write( io, str::Char ) = plug_write( io, str, plug='\t' )

end_write{S <: AbstractString}( io, str::S ) = plug_write( io, str, plug='\n' )
end_write( io, str::Char ) = plug_write( io, str, plug='\n' )

function writepsi( io, psi::Bio.Intervals.Interval{PsiEntry}, meta::AbstractVector, coord )
   for s in meta
      tab_write( io, s )
   end
   tab_write( io, psi.metadata.gene )
   tab_write( io, string(psi.metadata.node) )
   tab_write( io, coord )
   tab_write( io, string(psi.metadata.psi) )
   tab_write( io, string(psi.metadata.ci_width) )
   tab_write( io, string(psi.metadata.total) )
   tab_write( io, psi.metadata.ci_range )
   tab_write( io, psi.metadata.complexity )
   tab_write( io, string(psi.metadata.entropy) )
   tab_write( io, psi.metadata.inc )
   end_write( io, psi.metadata.exc )
end

main()
