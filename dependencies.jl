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
map( check_and_install, adds )
check_and_install("SuffixArrays", checkout=true)
check_and_install("IndexableBitVectors", checkout=true)

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
