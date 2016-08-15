#!/usr/bin/env julia

function check_and_install( pkg; clone=false, checkout=false )
   print( STDERR, "Checking $pkg ... " )
   pkgname = clone ? basename(pkg) |> x->split(x, ".jl.git", keep=false)[1] : pkg
   ver = Pkg.installed(pkgname)
   if ver != nothing
      println( STDERR, "Found version $ver" )
   else
      println( STDERR, "Trying to install $pkg ..." )
      if clone
         Pkg.clone(pkg)
      else
         Pkg.add(pkg)
      end
      if checkout
         Pkg.checkout(pkg)
      end
   end
end

adds = [ "DataStructures",
         "ArgParse",
         "FMIndexes", 
         "IntArrays",
         "IntervalTrees",
         "Libz",
         "StatsBase",
         "Distributions",
         "Glob",
         "Bio",
         "BaseTestNext" ]

Pkg.update()
#check_and_install("https://github.com/BioJulia/IndexableBitVectors.jl.git", clone=true, checkout=true )
#check_and_install("https://github.com/quinnj/SuffixArrays.jl.git", clone=true, checkout=true )
check_and_install("SuffixArrays", checkout=true)
check_and_install("IndexableBitVectors", checkout=true)
map( check_and_install, adds )


using DataStructures
using ArgParse
using Bio
using FMIndexes
using IntArrays
using IntervalTrees
using BufferedStreams
using Libz
using StatsBase
using Distributions
using Glob
using SuffixArrays
