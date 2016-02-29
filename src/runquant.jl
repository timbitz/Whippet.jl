#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

#using Bio.Seq
#using FMIndexes
#using Libz
using ArgParse
#using BufferedStreams
 
#include("types.jl")
#include("bio_nuc_safepatch.jl")
#include("refflat.jl")
#include("graph.jl")
#include("edges.jl")
#include("index.jl")
#include("align.jl")
#include("quant.jl")
#include("reads.jl")

dir = splitdir(@__FILE__)[1]

push!( LOAD_PATH, dir )
import SpliceGraphs
@everywhere using SpliceGraphs

function parse_cmd()
  s = ArgParseSettings(version="0.0.1", add_version=true)
  # TODO finish options...
  @add_arg_table s begin
    "filename.fastq[.gz]"
      arg_type = ASCIIString
      required = true
    "paired_mate.fastq[.gz]"
      arg_type = ASCIIString
      required = false
    "--no-tpm"
      help     = "Should tpm file be sent to output/prefix.tpm.gz? (default on)"
      action   = :store_true
    "--seed_len", "-K"
      help = "Seed length (default 18)"
      arg_type = Int
      default  = 18
    "--seed_try", "-M"
      help = "Number of failed seeds to try before giving up (default 3)"
      arg_type = Int
      default  = 3
    "--index", "-x"
      help = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = ASCIIString
      default = "$(dir)/../index/graph"
    "--out", "-o"
      help = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = ASCIIString
      default  = "$(dir)/../output"
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, "Loading splice graph index... $( args["index"] ).jls")
   @time const lib = open(deserialize, "$( args["index"] ).jls")

   println(STDERR, "Loading annotation index... $( args["index"] )_anno.jls")
   @time const anno = open(deserialize, "$( args["index"] )_anno.jls")

   const param = AlignParam() # defaults for now
   const quant = GraphLibQuant( lib, anno )
   const multi = Vector{Multimap}()

   parser = make_fqparser( fixpath(args["filename.fastq[.gz]"]) )

   if nprocs() > 1
      #include("align_parallel.jl")
      # Load Fastq files in chunks
      # Parallel reduction loop through fastq chunks 
      return #TODO
   else
      println(STDERR, "Processing reads...")
      @time mapped,total,readlen = process_reads!( parser, param, lib, quant, multi )
      readlen = round(Int, readlen)
      println(STDERR, "Finished $mapped mapped reads of length $readlen out of a total $total reads...")
   end

   # TPM_EM
   println(STDERR, "Calculating expression values...")
   calculate_tpm!( quant, readlen=readlen )
   @time iter = rec_gene_em!( quant, multi, sig=6, readlen=readlen )
   println(STDERR, "Finished calculating transcripts per million (TpM) after $iter iterations of EM...")

   if !args["no-tpm"]
      output_tpm( args["out"] * ".tpm.gz", lib, quant )
   end

   println(STDERR, "Assigning multi-mapping reads based on maximum likelihood estimate..")
   # Now assign multi to edges.
   @time assign_ambig!( quant, multi )

   println(STDERR, "Calculating effective lengths...")
   @time effective_lengths!( lib, quant, readlen - 19, min(readlen - param.score_min, 9-1) )
   @time bias_ave,bias_var = global_bias( quant )
   println(STDERR, "Global bias is $bias_ave +/- $bias_var ")
   println(STDERR, "Calculating maximum likelihood estimate of events..." )
   @time process_events( args["out"], lib, anno, quant )
   println(STDERR, "Whippet done." )
end

main()
