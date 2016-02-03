#!/usr/bin/env julia

function check_and_install( pkg )
   try
      print( STDERR, "Checking $pkg ... " )
      ver = Pkg.installed(pkg)
      println( STDERR, "Found version $ver" )
   catch
      println( STDERR, "Trying to install $pkg ..." )
      Pkg.add(pkg)
      Pkg.test(pkg)
   end
end

pkgs = [ "ArgParse", 
         "Bio", 
         "SuffixArrays", 
         "FMIndexes", 
         "IntArrays", 
         "IntervalTrees", 
         "GZip" ]
Pkg.update()
map( check_and_install, pkgs )

