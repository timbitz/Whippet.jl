#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Bio.Seq
using FMIndexes
using IntArrays
using GZip

include("kmerhash.jl")

K = 32
N = 100000000 # should be ~ 3billion kmers in the human genome

println(STDERR, "Allocating k-mer hash...")
gfmind = FMIndex[]

# replace N with A
function twobit_enc(seq)
    len = length(seq)
    ret = IntVector{2,UInt8}(len)
    for i in 1:len
        if seq[i] == DNA_N
            ret[i] = convert(UInt8, DNA_A)
        else
            ret[i] = convert(UInt8, seq[i])
        end
    end
    ret
end

function buildindex!( fhIter, indxArr; verbose=false )
    for r in fhIter
        immutable!(r.seq)
        println( STDERR, r )
        @time fm = FMIndex(twobit_enc(r.seq), 4, r=6)
        push!(indxArr, fm)
    end
    println( STDERR, "Finished building Index..." )
end 

# iterate through files in cd
for f in readdir()
   re = match(r"(\S+).(fa|fasta)(.gz)?", f)
  
   if re != nothing
      mainname = re.captures[1]
      if re.captures[3] != nothing #then gzipped
         println(STDERR, "Decompressing and Indexing $mainname...")
         #to_open = pipeline(`cat $f`, `gzip -dc`)
         continue 
      else
         println(STDERR, "Indexing $mainname...")
         to_open = f
      end
      # iterate through fasta entries
      @time buildindex!(open( to_open, FASTA ), gfmind)
   end
end


println()
println(STDERR, "Saving genome index...")
open(string("genome", ".jls"), "w+") do io
   @time serialize( io, gfmind )
end
gfmind = 0
gc()

println(STDERR, "Loading genome index...")
@time find = open(deserialize, "genome.jls")
