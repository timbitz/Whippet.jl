#__precompile__()

module SpliceGraphs

   import BufferedStreams
   import Bio
   import Bio.Seq
   import FMIndexes
   import IntervalTrees
   import Libz

   using BufferedStreams
   using Bio.Seq
   using FMIndexes
   using Libz
   using Base.Threads

   include("types.jl")
   include("bio_nuc_safepatch.jl")
   include("refflat.jl")
   include("graph.jl")
   include("edges.jl") 
   include("index.jl")
   include("align.jl")
   include("quant.jl")
   include("reads.jl")

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
          fasta_to_index,
          process_reads!,
          calculate_tpm!,
          rec_gene_em!,
          fixpath 
 
end
