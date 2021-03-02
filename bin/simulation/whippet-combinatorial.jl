#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/../VERSION")))

start = time_ns()
println( stderr, "Whippet $ver loading and compiling... " )

Pkg.activate(dir * "/..")

using Libz
using ArgParse
using StatsBase
using BioSequences
using Combinatorics

using Whippet

function parse_cmd()
  s = ArgParseSettings()
  # TODO finish options...
  @add_arg_table s begin
    "--index", "-x"
      help = "Prefix or full-name of index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = String
      default  = fixpath( "$(dir)/../../index/graph" )
    "--out", "-o"
      help = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = String
      default  = fixpath( "$(dir)/../../simul" )
    "--num-nodes", "-n"
      help = "Maximum number of consecutive nodes to output sliding combinatorial paths through"
      arg_type = Int64
      default  = 8
    "--num-genes", "-g"
      help = "Choose the first -n genes from the index to simulate (starting at offset -i)."
      arg_type = Int64
      default  = 100000
    "--offset", "-i"
      help = "Offset of gene number to start simulating from"
      arg_type = Int64
      default  = 1
    "--theoretical"
      action = :store_true
    "--verbose"
      help = "Print progress messages.."
      action = :store_true
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )

   println(stderr, "Loading splice graph index... $( args["index"] ).jls")
   @timer const lib = open(deserialize, "$( args["index"] ).jls")

   println(stderr, "Simulating combinatorial transcripts..")
   @timer simulate_genes( lib, output=args["out"],
                               gene_num=args["num-genes"],
                               offset=args["offset"],
                               node_num=args["num-nodes"],
                               verbose=args["verbose"],
                               theoretical=args["theoretical"] )
end

mutable struct SimulTranscript
   seq::Whippet.SGSequence
   nodes::Vector{UInt}
end

struct SimulGene
   trans::Vector{SimulTranscript}
   gene::String
end

istxstart( edge::EdgeType ) = edge == EDGETYPE_SL || edge == EDGETYPE_LS ? true : false
istxstop(  edge::EdgeType ) = edge == EDGETYPE_SR || edge == EDGETYPE_RS ? true : false

function collect_txstarts( sg::SpliceGraph )
   starts = Int[]
   for i in 1:length(sg.nodelen)
      if istxstart( sg.edgetype[i] )
         push!( starts, i )
      end
   end
   starts
end

function collect_txstops( sg::SpliceGraph )
   stops = Int[]
   for i in 2:length(sg.edgetype)
      if istxstop( sg.edgetype[i] )
         push!( stops, i-1 )
      end
   end
   stops
end

function combinatorial_paths( sg::SpliceGraph; max_nodes=10, verbose=false, theoretical=false )
   verbose && println("collecting starts & stops")
   starts = collect_txstarts( sg )
   stops  = collect_txstops( sg )
   paths  = Vector{Vector{Int16}}()
   verbose && println("going through $(length(starts)) starts and $(length(stops)) stops")
   isocnt = 0
   junc = Dict{Tuple{Int,Int},Int}()

   function add_juncs!( dict, path )
      for i in 1:length(path)-1
         dict[(path[i],path[i+1])] = 1
      end
   end

   function add_combinations!( paths, i, j, sub_range, val )
       for m in combinations( sub_range, val )
          for iv in reverse(i)
             pushfirst!( m, iv )
          end
          for jv in j
             push!( m, jv )
          end
          if isvalid_path( m, sg )
             isocnt += 1
             if theoretical
                add_juncs!( junc, m )
             else
                push!( paths, map(x->convert(Int16,x), m) )
             end
          end
       end
   end

   for i in starts, j in stops
      if i == j
         push!( paths, Int16[i] )
      elseif i == j - 1
         push!( paths, Int16[i,j] )
      elseif i < j
         range = i+1:j-1
         for k in 0:length(range)
            verbose && println(" going through combinations of $range, $k")
            if length(range) >= max_nodes
               sub_range = first(range):(first(range)+(length(range)-max_nodes))
               for r in sub_range
                  verbose && println(" adding sub_combinations of $(r:r+max_nodes-1) with $(i:r-1) and $(r+max_nodes:j) flanking")
                  add_combinations!( paths, i:r-1, r+max_nodes:j, r:r+max_nodes-1, k )
               end
            else
               add_combinations!( paths, [i], [j], range, k )
            end
         end
      end
   end
   paths, isocnt, length(values(junc))
end

function isvalid_path( path::Vector{Int}, sg::SpliceGraph )
   for i in 1:length(path)-1
      if path[i]+1 == path[i+1] # adjacent nodes
         if sg.edgetype[path[i+1]] in (EDGETYPE_LS, EDGETYPE_SR)
            return false
         end
      else
         if !(sg.edgetype[path[i]+1] in (EDGETYPE_LS, EDGETYPE_LR, EDGETYPE_LL) &&
              sg.edgetype[path[i+1]] in (EDGETYPE_SR, EDGETYPE_LR, EDGETYPE_RR))
            return false
         end
      end
   end
   true
end


function simulate_genes( lib; output="simul_genes", gene_num=length(lib.graphs),
                              offset=1, node_num=10, verbose=false, theoretical=false )
   fastaout = open( output * ".fa.gz", "w" )
   gtfout   = open( output * ".gtf.gz", "w" )
   fastastr = ZlibDeflateOutputStream( fastaout )
   gtfstr = ZlibDeflateOutputStream( gtfout )

   isocnt  = 0
   junccnt = 0

   for g in offset:min(offset+gene_num-1, length(lib.graphs)) #sample( 1:length(lib.graphs), min( length(lib.graphs), gene_num ), replace=false )
      sgene = SimulGene( Vector{SimulTranscript}(), lib.names[g] )
      iso,junc = simulate_transcripts( fastastr, gtfstr, sgene, lib.graphs[g], lib.info[g],
                                             node_num=node_num, verbose=verbose, theoretical=theoretical )
      #output_nodes( nodesstr, lib.names[g], lib.info[g], lib.graphs[g] )
      isocnt  += iso
      junccnt += junc
   end
   println(stderr, "$isocnt Unique Isoforms Simulated...")
   println(stderr, "$junccnt Unique Exon-exon Junctions..")
   close( fastastr )
   close( fastaout )
   close( gtfstr )
   close( gtfout )
end

function simulate_transcripts( fastream, gtfstream, simul::SimulGene, sg::SpliceGraph, info::GeneInfo;
                               node_num=10, verbose=false, theoretical=false )
   verbose && println("simulating paths for length $(length(sg.nodelen))")
   paths,isocnt,junccnt = combinatorial_paths( sg, max_nodes=node_num, verbose=verbose, theoretical=theoretical )
   verbose && println("length $(length(paths))")
   theoretical && return isocnt,junccnt
   for s in paths
      txid = simul.gene * "_" * join( s, "-" )
      output_transcript( fastream, s, sg, header=txid )
      output_gtf( gtfstream, s, sg, info, gene_id=simul.gene, transcript_id=txid )
   end
   isocnt, junccnt
end

function output_transcript( stream, path, sg::SpliceGraph; header::String="PREFIX" )
   write( stream, '>' )
   write( stream, header )
   write( stream, "\n" )
   for n in path
      noderange = sg.nodeoffset[n]:(sg.nodeoffset[n]+sg.nodelen[n]-1)
      for i in sg.seq[ noderange ]
         write( stream, convert(Char, i) )
      end
   end
   write( stream, "\n" )
end

function output_gtf( stream, path, sg::SpliceGraph, info::GeneInfo; gene_id::String="GENE", transcript_id::String="TXID" )
   start,stop = 0,0
   i = 1
   while i <= length(path)
      if i < length(path) && (sg.nodecoord[path[i+1]] == sg.nodecoord[path[i]]+sg.nodelen[path[i]] ||
                              sg.nodecoord[path[i+1]]+sg.nodelen[path[i+1]] == sg.nodecoord[path[i]])
         if info.strand
            start = start > 0 ? start : sg.nodecoord[path[i]]
            stop  = sg.nodecoord[path[i+1]]+sg.nodelen[path[i+1]]-1
         else
            start = sg.nodecoord[path[i+1]]
            stop  = stop > 0 ? stop : sg.nodecoord[path[i]]+sg.nodelen[path[i]]-1
         end
      else
         start = start > 0 ? start : Int(sg.nodecoord[path[i]])
         stop  = stop  > 0 ? stop  : Int(sg.nodecoord[path[i]]+sg.nodelen[path[i]]-1)

         write( stream, info.name * "\tWHIPPET_COMBINATORIAL\texon\t" )
         write( stream, string(start) * "\t" )
         write( stream, string(stop) * "\t\.\t" )
         write( stream, info.strand ? '+' : '-' )
         write( stream, "\t\.\tgene_id \"" * gene_id * "\"\; " )
         write( stream, "transcript_id \"" * transcript_id * "\"\;\n" )

         start,stop = 0,0
      end
      i += 1
   end
end

function output_nodes( stream, gname::String, info, sg::SpliceGraph)
   for n in 1:length(sg.nodelen)
      coord = "$(info.name):$(sg.nodecoord[n])-$(sg.nodecoord[n] + sg.nodelen[n] - 1)"
      write( stream, "$gname\t$n\t$coord\t$(info.gene)\n" )
   end
end

@timer main()
