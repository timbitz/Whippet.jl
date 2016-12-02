#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

tic()
println( STDERR, "Whippet $ver loading and compiling... " )

using ArgParse

push!( LOAD_PATH, dir * "/../src" )
using SpliceGraphs
using Libz

#=
using DataStructures
using BufferedStreams
using Bio.Seq
using FMIndexes
using IntArrays
using IntervalTrees
using Libz
using Distributions
using Requests

include("$dir/../src/types.jl")
include("$dir/../src/timer.jl")
include("$dir/../src/sgsequence.jl")
include("$dir/../src/fmindex_patch.jl")
include("$dir/../src/refset.jl")
include("$dir/../src/graph.jl")
include("$dir/../src/edges.jl")
include("$dir/../src/index.jl")
include("$dir/../src/align.jl")
include("$dir/../src/quant.jl")
include("$dir/../src/reads.jl")
include("$dir/../src/ebi.jl")
include("$dir/../src/paired.jl")
include("$dir/../src/events.jl")
include("$dir/../src/io.jl")
include("$dir/../src/diff.jl")
=#

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--kmer", "-k"
      help = "Kmer size to use for exon-exon junctions (default 9)"
      arg_type = Int64
      default  = 9
    "--fasta"
      help = "File containg the genome in fasta, one entry per chromosome [.gz]"
      arg_type = String
      required = true
    "--flat"
      help = "Gene annotation file in RefFlat format"
      arg_type = String
    "--gtf"
      help = "Gene anotation file in GTF format"
      arg_type = String
    "--index"
      help = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = String
      default = "$dir/../index/graph"
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()
   
   println(STDERR, " $( round( toq(), 6 ) ) seconds." )

   if args["gtf"] == nothing && args["flat"] == nothing
      println(STDERR, "ERROR: Must supply gene annotation file using `--gtf` or `--flat`!!")
      exit(1)
   elseif args["gtf"] != nothing
      annotype = "gtf"
      annotxt  = "GTF"
   else
      annotype = "flat"
      annotxt  = "Refflat"
   end

   println(STDERR, "Loading $annotxt file...")
   flat = fixpath( args[annotype] )
   fh = open( flat , "r")
   if isgzipped( flat )
      fh = fh |> x->ZlibInflateInputStream(x, reset_on_end=true)
   end
   @timer ref = annotype == "gtf" ? load_gtf(fh) : load_refflat(fh)

   println(STDERR, "Indexing transcriptome...")
   @timer graphome = fasta_to_index( fixpath( args["fasta"] ), ref, kmer=args["kmer"] )

   #=
   println(STDERR, "Saving Annotations...")
   open("$(args["index"])_anno.jls", "w") do fh
      @timer serialize(fh, ref)
   end
=#
   println(STDERR, "Serializing splice graph index...")
   open("$(args["index"]).jls", "w") do io
      @timer serialize( io, graphome )
   end

   println(STDERR, "Saving splice graph index...")
   @timer JLD.save("$(args["index"]).jld", "lib", graphome)
   
#=   indexpath = fixpath( args["index"] )
   println(STDERR, "Loading splice graph index... $( indexpath ).jls")
   @timer const lib = open(deserialize, "$( indexpath ).jls")

   println(STDERR, "Loading annotation index... $( indexpath )_anno.jls")
   @timer const anno = open(deserialize, "$( indexpath )_anno.jls")

   @timer res = JLD.load("$(args["index"]).jld")
   println(STDERR, "Whippet $ver done.")=#
end

main()
