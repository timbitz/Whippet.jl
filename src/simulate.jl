#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using ArgParse

dir = splitdir(@__FILE__)[1]

push!( LOAD_PATH, dir * "/../src" )
import SpliceGraphs
using SpliceGraphs

function parse_cmd()
  s = ArgParseSettings(version="Whippet v0.0.1-dev", add_version=true)
  # TODO finish options...
  @add_arg_table s begin
    "--index", "-x"
      help = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = ASCIIString
      default  = fixpath( "$(dir)/../index/graph" )
    "--out", "-o"
      help = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = ASCIIString
      default  = fixpath( "$(dir)/../output" )
    "--max-complexity", "-c"
      help = "What is the maximum complexity we should allow?"
      arg_type = Int64
      default  = 8
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, "Loading splice graph index... $( args["index"] ).jls")
   @time const lib = open(deserialize, "$( args["index"] ).jls")

   println(STDERR, "Loading annotation index... $( args["index"] )_anno.jls")
   @time const anno = open(deserialize, "$( args["index"] )_anno.jls")

   println(STDERR, "Simulating alternative transcripts of $( args["max-complexity"] ) max-complexity..")
   @time simulate_genes( lib, anno, args["max-complexity"] )

end

type SimulTranscript
   seq::SGSequence
   nodes::Vector{UInt}
end

immutable SimulGene
   trans::Vector{SimulTranscript}
   gene::ASCIIString
   complexity::Int
end

istxstart_left( edge::EdgeType ) = edge == EDGETYPE_SL || edge == EDGETYPE_LS ? true : false
istxstop_left(  edge::EdgeType ) = edge == EDGETYPE_SR || edge == EDGETYPE_RS ? true : false

function collect_nodes!( st::SimulTranscript, sg::SpliceGraph, r::UnitRange; skip=Vector{Int}, fromdonor=false, isdependent=false )
   n = r.start
   while n <= r.stop
      if     fromdonor && n in skip
         isdependent = true
      elseif fromdonor && (istxstart( sg.edgetype[n] ) || istxstop( sg.edgetype[n+1] ))
         isdependent = true
      elseif isdependent && !( sg.edgetype[n] in (EDGETYPE_RR, EDGETYPE_LR, EDGETYPE_SR) )
      else
         noderange = sg.nodeoffset[n]:(sg.nodeoffset[n]+sg.nodelen[n]-1)
         st.seq *= sg.seq[ noderange ]
         push!( st.nodes, n )
         fromdonor = sg.edgetype[n+1] in (EDGETYPE_LS, EDGETYPE_LR, EDGETYPE_LL) ? true : false
         isdependent = false
      end
      n += 1
   end
end

function simulate_genes( lib, anno, max_comp, output="simul_genes" )
   for g in 1:length(lib.graphs)
      comp = min( max_comp, length(lib.graphs[g].nodelen) - 2 )
      sgene = SimulGene( Vector{SimulTranscript}(), lib.names[g], rand(1:comp) )
      simulate_transcripts!( sgene, lib.graphs[g] )
      output_transcripts( fastaout, sgene, lib.graphs[g] )
      output_nodes( nodesout, lib.graphs[g] )
   end
end

function simulate_transcripst!( simul::SimulGene, sg::SpliceGraph )
   trans = SimulTranscript( sg"", Vector{UInt}() )
   collect_nodes!( trans, sg, 1:length(sg.nodelen) )
   push!( simul.trans, trans )
   if length(sg.nodelen) > 2
      mod_start = rand( 1:(length(sg.nodelen) - simul.complexity) )
      for i in 1:simul.complexity
         combset = collect( combinations( mod_start:(mod_start+simul.complexity), i ) )
         for s in combset
            trans = SimulTranscript( sg"", Vector{UInt}() )
            collect_nodes!( trans, sg, 1:length(sg.nodelen), skip=s )
            push!( simul.trans, trans )
         end
      end
   end
end

function output_transcripts( simul::SimulGene, sg::SpliceGraph )
   for t in simul.trans
      
   end
end

main()
