__precompile__()

module SpliceGraphs

   import DataStructures
   import BufferedStreams
   import Bio
   import Bio.Seq
   import FMIndexes
   import IntArrays
   import IntervalTrees
   import Libz

   using DataStructures
   using BufferedStreams
   using Bio.Seq
   using FMIndexes
   using IntArrays
   using IntervalTrees
   using Libz

   include("types.jl")
   include("bio_nuc_safepatch.jl")
   include("refflat.jl")
   include("graph.jl")
   include("edges.jl") 
   include("index.jl")
   include("align.jl")
   include("quant.jl")
   include("reads.jl")
   include("events.jl")

   if VERSION >= v"0.5.0-dev"
      using Base.Threads
      include("threaded.jl")
   end

   export SpliceGraph,
          SpliceGraphQuant,
          GraphLib,
          GraphLibQuant,
          SGNode,
          Edges,
          EdgeType,
          Kmer,
          SGKmer,
          SGNucleotide,
          SGNucleotideSequence,
          NucleotideSequence,
          Refset, Reftx, Refgene,
          Multimap,
          AlignParam,
          ungapped_align,
          make_fqparser,
          isgzipped,
          load_refflat,
          fasta_to_index,
          process_reads!,
          calculate_tpm!,
          rec_gene_em!,
          assign_ambig!,
          fixpath,
          effective_lengths!,
          global_bias,
          process_events 
end
