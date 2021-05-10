#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = readline(open(dir * "/../VERSION"))

start = time_ns()
println( stderr, "Whippet $ver loading and compiling... " )

Pkg.activate(dir * "/..")

using ArgParse
using Libz

using Whippet

function parse_cmd()
   s = ArgParseSettings()
   @add_arg_table s begin
    "--tpm"
      help     = "File containing two columns, isoform names, and transcripts-per-million (TpM)"
      arg_type = String
    "--index", "-x"
      help     = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph[.jls])"
      arg_type = String
      default  = fixpath( "$dir/../../index/graph" )
    "--out", "-o"
      help     = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = String
      default  = fixpath( "./output" )
   end
   return parse_args(s)
end

function main()

   args = parse_cmd()

   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )

   indexpath = fixpath( args["index"] )
   println(stderr, "Loading splice graph index... $( indexpath ).jls")
   @timer const lib = open(deserialize, "$( indexpath ).jls")

   data = load_tpm_file( fixpath(args["tpm"]) )

   tpms = Vector{Vector{Float64}}(length(lib.graphs))
   allocate_tpms!( tpms, lib )
   fill_tpms!( tpms, lib, data )

   io = open( args["out"], "w" )
   stream = ZlibDeflateOutputStream( io )

   output_tpms( stream, tpms, lib )

   close(stream)
   close(io)

   println(stderr, "Whippet $ver done." )
end

function load_tpm_file( filename::String, header=true )
   tpmfile = Dict{String,Float64}()
   fh = open( filename )
   for l in eachline(fh)
      if header
         header = false
         continue
      end
      s = split(l)
      tpmfile[String(s[1])] = parse(Float64, s[2])
   end
   tpmfile
end

function allocate_tpms!( tpms, lib )
   for i in 1:length(lib.graphs)
      tpms[i] = zeros(length(lib.graphs[i].annoname))
   end
   tpms
end

function fill_tpms!( tpms, lib, data )
   for i in 1:length(lib.graphs)
      for j in 1:length(lib.graphs[i].annoname)
         if haskey( data, lib.graphs[i].annoname[j] )
            tpms[i][j] = data[ lib.graphs[i].annoname[j] ]
         end
      end
   end
   tpms
end

function flanking_spliced( edges, node, name )
   upstream   = false
   downstream = false
   for i in node:-1:1
      if edges[i] in (EDGETYPE_RR, EDGETYPE_LR, EDGETYPE_SR)
         upstream = true
         #print("$name $(edges[i]):$i for n:$node")
         break
      elseif edges[i] in (EDGETYPE_LS, EDGETYPE_SL, EDGETYPE_RS)
         return false
      end
   end
   for i in node+1:length(edges)
      if edges[i] in (EDGETYPE_LL, EDGETYPE_LR, EDGETYPE_LS)
         downstream = true
         #println(" to $(edges[i]):$i")
         break
      elseif edges[i] in (EDGETYPE_SR, EDGETYPE_RS, EDGETYPE_SL)
         return false
      end
   end
   return (upstream && downstream)
end

function output_tpms( stream, tpms, lib )
   for i in 1:length(lib.graphs)
      for n in 1:length(lib.graphs[i].nodelen)
         if lib.graphs[i].edgetype[n]   in (EDGETYPE_RR, EDGETYPE_LR, EDGETYPE_LL) &&
            lib.graphs[i].edgetype[n+1] in (EDGETYPE_LL, EDGETYPE_LR, EDGETYPE_RR) &&
            flanking_spliced(lib.graphs[i].edgetype, n, lib.info[i].gene)

             tpm_total = 0.0
             tpm_node  = 0.0
             for j in 1:length(lib.graphs[i].annopath)
                curpath = lib.graphs[i].annopath[j]
                if n in first(curpath):last(curpath) # is within bounds
                   tpm_total += tpms[i][j]
                   if n in curpath
                      tpm_node += tpms[i][j]
                   end
                end
             end
             if tpm_total > 0.0
                psi = tpm_node / tpm_total
                write( stream, lib.info[i].gene * "\t" )
                write( stream, lib.info[i].name * ":" * string(lib.graphs[i].nodecoord[n]) *
                                                  "-" * string(lib.graphs[i].nodecoord[n]+lib.graphs[i].nodelen[n]-1) * "\t" )
                write( stream, string(n) * "\t" )
                write( stream, string(psi) * "\n" )
             end
         end
      end
   end
end

@timer main()
