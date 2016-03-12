#!/usr/bin/env julia

function check_and_install( pkg )
   print( STDERR, "Checking $pkg ... " )
   ver = Pkg.installed(pkg)
   if ver != nothing
      println( STDERR, "Found version $ver" )
   else
      println( STDERR, "Trying to install $pkg ..." )
      Pkg.add(pkg)
      #Pkg.test(pkg)
   end
end

pkgs = [ "DataStructures",
         "ArgParse", 
         "Bio", 
         "SuffixArrays", 
         "FMIndexes", 
         "IntArrays",
         "IntervalTrees",
         "BufferedStreams", 
         "Libz",
         "StatsBase" ]
Pkg.update()
map( check_and_install, pkgs )

