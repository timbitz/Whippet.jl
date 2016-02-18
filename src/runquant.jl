#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Bio.Seq
using FMIndexes

using ArgParse
 
include("types.jl")
include("bio_nuc_safepatch.jl")
include("refflat.jl")
include("graph.jl")
include("edges.jl")
include("index.jl")
include("align.jl")
include("quant.jl")
include("reads.jl")

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

   parser = make_fqparser( "../test/sample_01.fq.gz" )

   if nprocs() > 1
      #include("align_parallel.jl")
      # Load Fastq files in chunks
      # Parallel reduction loop through fastq chunks 
      return #TODO
   else
      println(STDERR, "Processing reads...")
      @time mapped,unmapped = process_reads!( parser, quant, multi, param, lib )
      println(STDERR, "Finished $mapped mapped reads out of a total $(mapped+unmapped) reads...")
   end

   # TPM_EM
   println(STDERR, "Calculating expression values...")
   calculate_tpm!( quant )
   @time iter = rec_gene_em!( quant, multi )
   println(STDERR, "Finished calculating transcripts per million (TpM) after $iter iterations of EM...")

   # Now assign multi to edges.

   # Iterate through Events and do PSI_EM 
   #calculate_psi( quant )
end

main()
