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
using BufferedStreams
using Libz


function parse_cmd()
  s = ArgParseSettings()
  @add_arg_table s begin
    "filename.psi.gz"
      arg_type = String
      required = true
    "--fasta"
      help     = "Genome sequence fasta file used with `--fasta` parameter when building the whippet index"
      arg_type = String
      required = true

  end
  return parse_args(s)
end

function remove_ambiguous!( seq, replace_with::DNA=DNA_A )
   for (pos, nuc) in each(isambiguous, seq)
      seq[pos] = replace_with
   end
end

function read_fasta( filename::String )
   genome = Dict{String,ReferenceSequence}()

   if isgzipped( filename )
      println(stderr, "Decompressing and Indexing $filename...")
      to_open = open( filename ) |> x->ZlibInflateInputStream(x, reset_on_end=true)
   else
      println(stderr, "Indexing $filename...")
      to_open = open( filename ) |> BufferedInputStream
   end
   reader = FASTA.Reader( to_open )

   for entry in reader
      println(stderr, "Processing $(identifier(entry))")
      seq = sequence(entry)
      remove_ambiguous!(seq)
      genome[identifier(entry)] = ReferenceSequence(seq)
   end

   genome
end

function parse_coord( coord::String )
   chr,inter = split(coord, ':')
   if '-' in inter
      left,right = map(x->parse(Int,x), split(inter, '-'))
      return chr,left,right
   else
      middle = parse(Int, inter)
      return chr,middle
   end
end

function main()
   args  = parse_cmd()
   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )

   genome = read_fasta( args["fasta"] )

   for l in eachline( ZLibInflateInputStream(open(args["filename.psi.gz"], "r")) )
      spl = split(l, '\t')
      ci, si, ti = 3,4,5 #coord, strand, type
      if (spl[ti] in ("XDI", "XDE") && spl[si] == '+') ||
         (spl[ti] in ("XAI", "XAE") && spl[si] == '-')

         name,coord = parse_coord(spl[ci])
         
      else
         name,left,right = parse_coord(spl[ci])
         spl[ti] 
      end

   end

   println(stderr, "Whippet $ver done." )
end

@timer main()
