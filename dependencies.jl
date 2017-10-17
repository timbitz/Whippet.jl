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
         "IntArrays",
         "FMIndexes",
         "Libz",
         "StatsBase",
         "Distributions",
         "Glob",
         "BioSymbols",
         "BioSequences" ]

checkouts = [ "IndexableBitVectors",
              "IntervalTrees" ]

tic()
Pkg.update()
map( x->check_and_install(x, checkout=true), checkouts )
map( check_and_install, adds )

println( STDERR, "INFO: Loading and precompiling dependencies... " )

println( STDERR, "INFO: using DataStructures" )
using DataStructures
println( STDERR, "INFO: using ArgParse" )
using ArgParse
println( STDERR, "INFO: using Bio" )
using BioSymbols
using BioSequences
println( STDERR, "INFO: using FMIndexes" )
using FMIndexes
println( STDERR, "INFO: using IntArrays" )
using IntArrays
println( STDERR, "INFO: using IntervalsTrees" )
using IntervalTrees
println( STDERR, "INFO: using BufferedStreams" )
using BufferedStreams
println( STDERR, "INFO: using Libz" )
using Libz
println( STDERR, "INFO: using StatsBase" )
using StatsBase
println( STDERR, "INFO: using Distributions" )
using Distributions
println( STDERR, "INFO: using Glob" )
using Glob
println( STDERR, "INFO: using SuffixArrays" )
using SuffixArrays
#println( STDERR, "INFO: using Requests" )
#using Requests

println( STDERR, "INFO: Loading and precompiling whippet... " )

const dir = abspath( splitdir(@__FILE__)[1] )
push!( LOAD_PATH, dir * "/src" )
using Whippet
toc()
