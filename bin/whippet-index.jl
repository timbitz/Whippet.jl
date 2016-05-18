#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

const ver = chomp(readline(open(splitdir(@__FILE__)[1] * "/VERSION")))
const dir = abspath( splitdir(@__FILE__)[1] * "/../src" )

tic()
println( STDERR, "Whippet $ver loading and compiling... " )

using DataStructures
using IntervalTrees
using Bio.Seq
using FMIndexes
using IntArrays
using Libz

using ArgParse

include("$dir/types.jl")
include("$dir/bio_nuc_safepatch.jl")
include("$dir/refflat.jl")
include("$dir/graph.jl")
include("$dir/edges.jl")
include("$dir/index.jl")
include("$dir/timer.jl")

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--kmer", "-k"
      help = "Kmer size to use for exon-exon junctions (default 9)"
      arg_type = Int64
      default  = 9
    "--fasta"
      help = "File containg the genome in fasta, one entry per chromosome [.gz]"
      arg_type = ASCIIString
      required = true
    "--flat"
      help = "Gene annotation file in RefFlat format"
      arg_type = ASCIIString
      required = true
    "--index"
      help = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = ASCIIString
      default = "$dir/../index/graph"
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()
   
   println(STDERR, " $( round( toq(), 6 ) ) seconds." )

   println(STDERR, "Loading Refflat file...")
   flat = fixpath( args["flat"] )
   fh = open( flat , "r")
   if isgzipped( flat )
      fh = fh |> x->ZlibInflateInputStream(x, reset_on_end=true)
   end
   @timer ref = load_refflat(fh)

   println(STDERR, "Indexing transcriptome...")
   @timer graphome = fasta_to_index( fixpath( args["fasta"] ) , ref, kmer=args["kmer"] )

   println(STDERR, "Saving Annotations...")
   open("$(args["index"])_anno.jls", "w+") do fh
      @timer serialize(fh, ref)
   end

   println(STDERR, "Saving splice graph index...")
   open("$(args["index"]).jls", "w+") do io
      @timer serialize( io, graphome )
   end

   println(STDERR, "Whippet $ver done.")
end

main()
