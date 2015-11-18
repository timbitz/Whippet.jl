#!/usr/bin/env julia

pkgs = ["ArgParse", "FMIndexes", "IntArrays", "IntervalTrees", "GZip"]
map( Pkg.add, pkgs )
Pkg.update()
map( Pkg.test, pkgs )

using FMIndexes
using IntArrays
using ArgParse
using IntervalTrees
using GZip
