#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Bio.Seq
using FMIndexes
using IntArrays

using ArgParse
 
include("types.jl")
include("index.jl")
include("bio_nuc_safepatch.jl")
include("refflat.jl")
include("graph.jl")
include("edges.jl")
include("align.jl")

function main()
   println(STDERR, "Loading Refflat file...")
   fh = open("$(pwd())/../genome/genes.flat", "r")
   @time ref = load_refflat(fh)
   close(fh)

   genomedir = string(pwd(), "/../genome")
   println(STDERR, "Indexing transcriptome...")
   graphome = fasta_to_index( genomedir, refflat )

   println(STDERR, "Saving Annotations...")
   open("$(pwd())/../index/flat.jls", "w+") do fh
      @time serialize(fh, ref)
   end

   println(STDERR, "Saving splice graph index...")
   open("$(pwd())/../index/graph.jls", "w+") do io
      @time serialize( io, graphome )
   end

end

main()
