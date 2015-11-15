#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Bio.Seq
using FMIndexes
using IntArrays
using GZip

include("kmerhash.jl")

include("index.jl")

K = 32
N = 100000000 # should be ~ 3billion kmers in the human genome

println(STDERR, "Allocating k-mer hash...")
gfmind = FMIndex[]
genome = nothing


# iterate through files in cd

function fasta_to_index( dir )
   genome = nothing
   for f in readdir(dir)
      re = match(r"(\S+).(fa|fasta)(.gz)?", f)
  
      if re != nothing
         mainname = re.captures[1]
         if re.captures[3] != nothing #then gzipped
            println(STDERR, "Decompressing and Indexing $mainname...")
            #to_open = pipeline(`cat $f`, `gzip -dc`)
            continue 
         else
            println(STDERR, "Indexing $mainname...")
            to_open = string(dir, "/$f")
         end
         # iterate through fasta entries
         #@time buildindex!(open( to_open, FASTA ), gfmind)
         genome = @time single_index!(open( to_open, FASTA ))
      end
   end
   genome
end

function main()
   genomedir = string(pwd(), "/genome")   
   println(STDERR, "Indexing genome...")
   genome = fasta_to_index( genomedir )

   println()
   println(STDERR, "Saving genome index...")
   open("$(pwd())/index/genome2.jls", "w+") do io
   #@time serialize( io, gfmind )
      @time serialize( io, genome )
   end

   println(STDERR, "Loading genome index...")
   @time find = open(deserialize, "$(pwd())/index/genome2.jls")
end

main()
