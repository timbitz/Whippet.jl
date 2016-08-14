#!/usr/bin/env julia

function check_and_install( pkg; clone=false )
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
   end
end

clones = [ "https://github.com/dcjones/Switch.jl.git",
           "https://github.com/BioJulia/IndexableBitVectors.jl.git",
           "https://github.com/quinnj/SuffixArrays.jl.git" ]


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
map( x->check_and_install(x,clone=true), clones )


using DataStructures
using Switch
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
