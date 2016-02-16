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

function parse_cmd()
  s = ArgParseSettings()
  # TODO finish options...
  @add_arg_table s begin
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

   println(STDERR, "Loading splice graph index...")
   @time lib = open(deserialize, "$(pwd())/index/graph.jls")

   println(STDERR, "Loading annotation index...")
   @time anno = open(deserialize, "$(pwd())/../index/graph_anno.jls")

   const ap = AlignParam() # defaults for now

   # Load Fastq files in chunks
   # Parallel reduction loop through fastq chunks 

   # Add up non-multi hit gene counts and splice graph quants
   # Make Multihit objects
   # TPM_EM
   # Iterate through Events and do PSI_EM 
end

main()
