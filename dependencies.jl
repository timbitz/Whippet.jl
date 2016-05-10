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
         "BufferedStreams",
         "Libz",
         "Bio",
         "StatsBase",
         "Distributions",
         "Glob" ]

Pkg.update()
map( x->check_and_install(x,clone=true), clones )
map( check_and_install, adds )

Pkg.pin("Bio", v"0.1")
Pkg.pin("BufferedStreams", v"0.0.3")

using DataStructures
using Switch
using ArgParse
using Bio
using SuffixArrays
using FMIndexes
using IntArrays
using IntervalTrees
using BufferedStreams
using Libz
using StatsBase
using Distributions
using Glob
