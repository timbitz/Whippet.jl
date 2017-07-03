#!/usr/bin/env julia

function check_and_install( pkg; clone=false, checkout=false )
   print( STDERR, "Checking $pkg ... " )
   pkgname = clone ? basename(pkg) |> x->split(x, ".jl.git", keep=false)[1] |> string : pkg
   try
      ver = Pkg.installed(pkgname)
      if !clone && ver == nothing
         error()
      end
      println( STDERR, "Found version $ver" )
   catch
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
         "SuffixArrays",
         "FMIndexes", 
         "IntArrays",
         "IntervalTrees",
         "Libz",
         "StatsBase",
         "Distributions",
         "Glob",
         "BioSymbols",
         "BioSequences",
         "Requests"]

tic()
Pkg.update()
check_and_install("IndexableBitVectors", checkout=true)
map( check_and_install, adds )

println( STDERR, "INFO: Loading and precompiling... " )

using DataStructures
using ArgParse
using BioSymbols
using BioSequences
using FMIndexes
using IntArrays
using IntervalTrees
using BufferedStreams
using Libz
using StatsBase
using Distributions
using Glob
using SuffixArrays
using Requests

const dir = abspath( splitdir(@__FILE__)[1] )
push!( LOAD_PATH, dir * "/src" )
using Whippet
toc()
