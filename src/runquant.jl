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

push!( LOAD_PATH, splitdir(@__FILE__)[1] )
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
    "--tpm"
      help     = "Should tpm file be sent to STDOUT? (default off)"
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
      default = "$(pwd())/../index/graph"
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, "Loading splice graph index...")
   @time const lib = open(deserialize, "$(pwd())/../index/graph.jls")

   println(STDERR, "Loading annotation index...")
   @time const anno = open(deserialize, "$(pwd())/../index/graph_anno.jls")

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

   if args["tpm"]
      for i in 1:length(lib.names)
         println(lib.names[i] * "\t" * string(quant.tpm[i]) ) 
      end
   end

   # Now assign multi to edges.
   assign_ambig!( quant, multi )

   @time effective_lengths!( lib, quant, readlen - 19, min(readlen - param.score_min, 9-1) )
   @time bias_ave,bias_var = global_bias( quant )
   println("Calculating global bias to $bias_ave +/- $bias_var ")
   @time process_events( "tmp.gz", lib, anno, quant )
end

main()
