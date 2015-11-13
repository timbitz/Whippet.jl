#!/usr/bin/env julia

pkgs = ["ArgParse", "Bio", "FMIndexes", "BitArrays", "GZip"]
map( x->Pkg.add(x), pkgs )
Pkg.update()
map( x->Pkg.test(x), pkgs )

using Bio.Seq
using FMIndexes
using BitArrays
using ArgParse
using GZip

println(STDERR, "looks okay!")
