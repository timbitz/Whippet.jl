#!/usr/bin/env julia

pkgs = ["ArgParse", "Bio", "FMIndexes", "IntArrays", "GZip"]
map( Pkg.add, pkgs )
Pkg.update()
map( Pkg.test, pkgs )

using Bio.Seq
using FMIndexes
using IntArrays
using ArgParse
using GZip

println(STDERR, "looks okay!")
