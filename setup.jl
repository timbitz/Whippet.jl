#!/usr/bin/env julia

pkgs = ["ArgParse", "Bio", "Glob", "GZip"]
map( x->Pkg.add(x), pkgs )
Pkg.update()

using Bio.Seq
using ArgParse
using Glob
using GZip
