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
include("reads.jl")
include("align.jl")
include("quant.jl")

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
   @time const lib = open(deserialize, "$(pwd())/../index/graph.jls")

   println(STDERR, "Loading annotation index...")
   @time const anno = open(deserialize, "$(pwd())/../index/graph_anno.jls")

   const param = AlignParam() # defaults for now
   const quant = GraphLibQuant( lib, anno )
   const multi = Vector{Multimap}()

   if nprocs() > 1
      include("align_parallel.jl")
      # Load Fastq files in chunks
      # Parallel reduction loop through fastq chunks 
   else
      parser = make_fqparser( "../test/sample_01.fq.gz" )
      reads  = allocate_chunk( parser, 10000 )
      while !done(parser)
         read_chunk!( reads, parser )
         for r in reads
            align = ungapped_align( param, lib, r )
            if !isnull( align )
               if length( get(align) ) > 1
                  push!( multi, Multimap( get(align) ) )
               else
                  count!( quant, get(align) )
               end
            end
         end
      end # end while
   end

   # TPM_EM
   # Iterate through Events and do PSI_EM 
end

main()
