__precompile__()

module Whippet

importall Bio

using Bio.Seq
using DataStructures
using BufferedStreams
using FMIndexes
using IntArrays
using IntervalTrees
using Libz
using Distributions
using Requests

include("types.jl")
include("timer.jl")
include("sgkmer.jl")
include("fmindex_patch.jl")
include("refset.jl")
include("graph.jl")
include("edges.jl") 
include("index.jl")
include("align.jl")
include("quant.jl")
include("reads.jl")
include("ebi.jl")
include("paired.jl")
include("events.jl")
include("io.jl")
include("diff.jl")

#=if VERSION >= v"0.5.0"
   using Base.Threads
   include("threaded.jl")
end=#

export SpliceGraph,
       SpliceGraphQuant,
       GraphLib,
       GraphLibQuant,
       GeneInfo,
       TxInfo,
       SGNode,
       Edges,
       EdgeType,
       EDGETYPE_LS,
       EDGETYPE_SL,
       EDGETYPE_RS,
       EDGETYPE_SR,
       EDGETYPE_LR,
       EDGETYPE_LL,
       EDGETYPE_RR,
       Kmer,
       SGKmer,
       SGNucleotide,
       SGNucleotideSequence,
       SGSequence,
       @sg_str,
       reverse_complement,
       kmer,
       NucleotideSequence,
       RefSet, RefTx, RefGene,
       Multimap,
       AlignParam,
       ungapped_align,
       make_fqparser,
       make_http_fqparser,
       isgzipped,
       EBIResponse,
       ident_to_fastq_url,
       load_refflat,
       load_gtf,
       fasta_to_index,
       process_reads!,
       process_paired_reads!,
       calculate_tpm!,
       rec_gene_em!,
       gene_em!,
       output_tpm,
       output_stats,
       output_junctions,
       assign_ambig!,
       effective_lengths!,
       global_bias,
       process_events,
       fixpath,
       @timer,
       timer_print,
       open_streams,
       process_psi_files
end
