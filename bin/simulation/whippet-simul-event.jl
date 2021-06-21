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

using Whippet

function parse_cmd()
  s = ArgParseSettings()
  # TODO finish options...
  @add_arg_table s begin
    "--index", "-x"
      help = "Prefix for index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = String
      default  = fixpath( "$(dir)/../../index/graph" )
    "--out", "-o"
      help = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = String
      default  = fixpath( "$(dir)/../../simul" )
    "--gtf"
      help = "Gene anotation file in GTF format"
      arg_type = String
    "--max-complexity", "-k"
      help = "Max complexity KMax allowed per gene (random K for each gene between K1:KMax)?"
      arg_type = Int64
      default  = 8
    "--num-genes", "-n"
      help = "Randomly sample at most -n number of genes from the index to simulate."
      arg_type = Int64
      default  = 5000
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )

   println(stderr, "Loading splice graph index... $( args["index"] ).jls")
   @timer const lib = open(deserialize, "$( args["index"] ).jls")

   println(stderr, "Simulating alternative transcripts of $( args["max-complexity"] ) max-complexity..")
   @timer simulate_genes( lib, anno, args["max-complexity"], output=args["out"], gene_num=args["num-genes"] )

end

type SimulTranscript
   seq::Whippet.SGSequence
   nodes::Vector{UInt}
end

immutable SimulGene
   trans::Vector{SimulTranscript}
   gene::String
   complexity::Int
end

istxstart( edge::EdgeType ) = edge == EDGETYPE_SL || edge == EDGETYPE_LS ? true : false
istxstop(  edge::EdgeType ) = edge == EDGETYPE_SR || edge == EDGETYPE_RS ? true : false

function collect_nodes!( st::SimulTranscript, sg::SpliceGraph, r::UnitRange; skip=Vector{Int}(), fromdonor=false, isdependent=false )
   n = r.start
   while n <= r.stop
      if     fromdonor && n in skip
         isdependent = true
      elseif fromdonor && istxstart( sg.edgetype[n] )
         isdependent = true
      elseif fromdonor && istxstop( sg.edgetype[n+1] ) && n < length(sg.nodelen)
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

function simulate_genes( lib, anno, max_comp; output="simul_genes", gene_num=length(lib.graphs) )
   fastaout = open( output * ".fa.gz", "w" )
   nodesout = open( output * ".node.gz", "w" )
   fastastr = ZlibDeflateOutputStream( fastaout )
   nodesstr = ZlibDeflateOutputStream( nodesout )

   for g in sample( 1:length(lib.graphs), min( length(lib.graphs), gene_num ), replace=false )
      comp = min( max_comp, length(lib.graphs[g].nodelen) - 2 )
      sgene = SimulGene( Vector{SimulTranscript}(), lib.names[g], rand(1:comp) )
      simulate_transcripts!( sgene, lib.graphs[g] )
      output_transcripts( fastastr, sgene, lib.graphs[g] )
      output_nodes( nodesstr, lib.names[g], anno.geneset[lib.names[g]].info, lib.graphs[g] )
   end
   close( fastastr )
   close( fastaout )
   close( nodesstr )
   close( nodesout )
end

function simulate_transcripts!( simul::SimulGene, sg::SpliceGraph )
   trans = SimulTranscript( dna"", Vector{UInt}() )
   collect_nodes!( trans, sg, 1:length(sg.nodelen) )
   push!( simul.trans, trans )
   if length(sg.nodelen) > 2
      mod_start = rand( 1:(length(sg.nodelen) - simul.complexity) )
      for i in 1:simul.complexity
         combset = collect( combinations( mod_start:(mod_start+simul.complexity-1), i ) )
         for s in combset
            trans = SimulTranscript( dna"", Vector{UInt}() )
            collect_nodes!( trans, sg, 1:length(sg.nodelen), skip=s )
            push!( simul.trans, trans )
         end
      end
   end
end

function output_transcripts( stream, simul::SimulGene, sg::SpliceGraph )
   prefix = simul.gene * "_K" * string(simul.complexity) * "_"
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
      coord = "$(info[1]):$(sg.nodecoord[n])-$(sg.nodecoord[n] + sg.nodelen[n] - 1)"
      write( stream, "$gname\t$n\t$coord\t$(info[2])\n" )
   end
end

main()
