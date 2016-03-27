#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using ArgParse
 
dir = splitdir(@__FILE__)[1]

push!( LOAD_PATH, dir * "/../src" )
import SpliceGraphs
using SpliceGraphs

function parse_cmd()
  s = ArgParseSettings(version="Whippet v0.0.1-dev", add_version=true)
  # TODO finish options...
  @add_arg_table s begin
    "filename.fastq[.gz]"
      arg_type = ASCIIString
      required = true
    "paired_mate.fastq[.gz]"
      arg_type = ASCIIString
      required = false
    "--index", "-x"
      help     = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = ASCIIString
      default  = fixpath( "$(dir)/../index/graph" )
    "--out", "-o"
      help     = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = ASCIIString
      default  = fixpath( "$(dir)/../output" )
    "--sam", "-s"
      help     = "Should SAM format be sent to stdout?"
      action   = :store_true
    "--seed_len", "-K"
      help     = "Seed length"
      arg_type = Int
      default  = 18
    "--seed_try", "-M"
      help     = "Number of failed seeds to try before giving up"
      arg_type = Int
      default  = 3
    "--junc-only", "-j"
      help     = "Only use junction reads, no internal exon reads will be considered."
      action   = :store_true
    "--no-tpm"
      help     = "Should tpm file be sent to output/prefix.tpm.gz? (default on)"
      action   = :store_true
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, "Loading splice graph index... $( args["index"] ).jls")
   @timer const lib = open(deserialize, "$( args["index"] ).jls")

   println(STDERR, "Loading annotation index... $( args["index"] )_anno.jls")
   @timer const anno = open(deserialize, "$( args["index"] )_anno.jls")

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
      @timer mapped,total,readlen = process_reads!( parser, param, lib, quant, multi, sam=args["sam"] )
      readlen = round(Int, readlen)
      println(STDERR, "Finished $mapped mapped reads of length $readlen out of a total $total reads...")
   end

   # TPM_EM
   println(STDERR, "Calculating expression values and MLE for $( length(multi) ) repetitive reads with EM...")
   calculate_tpm!( quant, readlen=readlen )
   @timer iter = rec_gene_em!( quant, multi, sig=6, readlen=readlen, max=1000 ) 
   println(STDERR, "Finished calculating transcripts per million (TpM) after $iter iterations of EM...")

   if !args["no-tpm"]
      output_tpm( args["out"] * ".tpm.gz", lib, quant )
   end

   println(STDERR, "Assigning multi-mapping reads based on maximum likelihood estimate..")
   # Now assign multi to edges.
   @timer assign_ambig!( quant, multi )

   println(STDERR, "Calculating effective lengths...")
   @timer effective_lengths!( lib, quant, readlen - 19, min(readlen - param.score_min, 9-1) )
   @timer bias_ave,bias_var = global_bias( quant )
   println(STDERR, "Global bias is $bias_ave +/- $bias_var ")
   println(STDERR, "Calculating maximum likelihood estimate of events..." )
   @timer process_events( args["out"] * ".psi.gz" , lib, quant, isnodeok=!args["junc-only"] )
   println(STDERR, "Whippet done." )
end

main()
