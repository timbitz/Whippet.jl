#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

tic()
println( STDERR, "Whippet $ver loading and compiling... " )

using Libz
using ArgParse
using StatsBase
using BioSequences
using Combinatorics

push!( LOAD_PATH, dir * "/../src" )
using Whippet

function parse_cmd()
  s = ArgParseSettings()
  # TODO finish options...
  @add_arg_table s begin
    "--index", "-x"
      help = "Prefix for index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = String
      default  = fixpath( "$(dir)/../index/graph" )
    "--out", "-o"
      help = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = String
      default  = fixpath( "$(dir)/../simul" )
    "--gtf"
      help = "Gene anotation file in GTF format"
      arg_type = String 
    "--num-genes", "-n"
      help = "Randomly sample at most -n number of genes from the index to simulate."
      arg_type = Int64
      default  = 5000
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, " $( round( toq(), 6 ) ) seconds" )

   println(STDERR, "Loading splice graph index... $( args["index"] ).jls")
   @timer const lib = open(deserialize, "$( args["index"] ).jls")

   println(STDERR, "Simulating combinatorial transcripts..")
   @timer simulate_genes( lib, output=args["out"], gene_num=args["num-genes"] )

end

type SimulTranscript
   seq::Whippet.SGSequence
   nodes::Vector{UInt}
end

immutable SimulGene
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

function combinatorial_paths( sg::SpliceGraph, max_nodes=20 )
   println("collecting starts & stops")
   starts = collect_txstarts( sg )
   stops  = collect_txstops( sg )
   paths  = Vector{IntSet}()
   println("going through $(length(starts)) starts and $(length(stops)) stops")

   function add_combinations!( paths, i, j, sub_range, val )
       for m in combinations( sub_range, val )
          for iv in reverse(i)
             unshift!( m, iv )
          end
          for jv in j
             push!( m, jv )
          end
          if isvalid_path( m, sg )
             push!( paths, IntSet(m) )
          end
       end
   end

   for i in starts, j in stops
      if i == j
         push!( paths, IntSet([i]) )
      elseif i == j - 1
         push!( paths, IntSet([i, j]))
      elseif i < j
         range = i+1:j-1
         for k in 1:length(range)
            println(" going through combinations of $range, $k")
            if length(range) >= max_nodes
               sub_range = first(range):(first(range)+(length(range)-max_nodes))
               for r in sub_range
                  println(" adding sub_combinations of $(r:r+max_nodes-1) with $(i:r-1) and $(r+max_nodes:j) flanking")
                  add_combinations!( paths, i:r-1, r+max_nodes:j, r:r+max_nodes-1, k )
               end
            else
               add_combinations!( paths, [i], [j], range, k )
            end
         end
      end
   end
   paths
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


function collect_nodes!( st::SimulTranscript, sg::SpliceGraph, path )
   for n in path
      noderange = sg.nodeoffset[n]:(sg.nodeoffset[n]+sg.nodelen[n]-1)
      st.seq *= sg.seq[ noderange ]
      push!( st.nodes, n )
   end
   st
end

function simulate_genes( lib; output="simul_genes", gene_num=length(lib.graphs) )
   fastaout = open( output * ".fa.gz", "w" ) 
   nodesout = open( output * ".node.gz", "w" )
   fastastr = ZlibDeflateOutputStream( fastaout )
   nodesstr = ZlibDeflateOutputStream( nodesout )
   
   for g in sample( 1:length(lib.graphs), min( length(lib.graphs), gene_num ), replace=false )
      sgene = SimulGene( Vector{SimulTranscript}(), lib.names[g] )
      simulate_transcripts!( sgene, lib.graphs[g] )
      output_transcripts( fastastr, sgene, lib.graphs[g] )
      output_nodes( nodesstr, lib.names[g], lib.info[g], lib.graphs[g] )
   end
   close( fastastr )
   close( fastaout )
   close( nodesstr )
   close( nodesout )
end

function simulate_transcripts!( simul::SimulGene, sg::SpliceGraph )
   println("simulating paths for length $(length(sg.nodelen))")
   paths = combinatorial_paths( sg )
   println("length $(length(paths))")
   for s in paths
      trans = SimulTranscript( dna"", Vector{UInt}() )
      collect_nodes!( trans, sg, s )
      push!( simul.trans, trans )
   end
end

function output_transcripts( stream, simul::SimulGene, sg::SpliceGraph )
   prefix = simul.gene * "_"
   used = Set()
   for t in simul.trans
      name = prefix * join( t.nodes, "-" )
      (name in used) && continue
      write( stream, ">$name\n" )
      for i in t.seq
         write( stream, Char(i) )
      end
      write( stream, "\n" )
      push!( used, name )
   end
end

function output_nodes( stream, gname::String, info, sg::SpliceGraph)
   for n in 1:length(sg.nodelen)
      coord = "$(info.name):$(sg.nodecoord[n])-$(sg.nodecoord[n] + sg.nodelen[n] - 1)"
      write( stream, "$gname\t$n\t$coord\t$(info.gene)\n" )
   end
end

main()
