#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Bio.Seq
using FMIndexes
using IntArrays
using Libz
using BufferedStreams

using ArgParse
 
include("types.jl")
include("bio_nuc_safepatch.jl")
include("refflat.jl")
include("graph.jl")
include("edges.jl")
include("index.jl")
#include("align.jl")

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--kmer", "-k"
      help = "Kmer size to use for exon-exon junctions (default 9)"
      arg_type = Int
      default  = 9
    "--fasta"
      help = "Directory or file containg the genome files in fasta"
      arg_type = ASCIIString
      required = true
    "--flat"
      help = "Gene annotation file in RefFlat format"
      arg_type = ASCIIString
      required = true
    "--index"
      help = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = ASCIIString
      default = "$(pwd())/../index/graph"
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, "Loading Refflat file...")
   fh = open( fixpath(args["flat"]) , "r")
   if isgzipped( args["flat"] )
      fh = fh |> ZlibInflateInputStream
   end
   @time ref = load_refflat(fh)

   genomedir = args["fasta"]
   println(STDERR, "Indexing transcriptome...")
   @time graphome = fasta_to_index( fixpath(genomedir) , ref )

   println(STDERR, "Saving Annotations...")
   open("$(args["index"])_anno.jls", "w+") do fh
      @time serialize(fh, ref)
   end

   println(STDERR, "Saving splice graph index...")
   open("$(args["index"]).jls", "w+") do io
      @time serialize( io, graphome )
   end

end

main()
