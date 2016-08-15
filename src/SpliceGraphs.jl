__precompile__()

module SpliceGraphs

import DataStructures
import BufferedStreams
import Bio
import FMIndexes
import IntArrays
import IntervalTrees
import Libz
import Distributions

using DataStructures
using BufferedStreams
using Bio.Seq
using FMIndexes
using IntArrays
using IntervalTrees
using Libz
using Distributions

include("types.jl")
include("timer.jl")
include("bio_nuc_safepatch.jl")
include("refset.jl")
include("graph.jl")
include("edges.jl") 
include("index.jl")
include("align.jl")
include("quant.jl")
include("reads.jl")
include("paired.jl")
include("events.jl")
include("io.jl")
include("diff.jl")

if VERSION >= v"0.5.0-dev"
   using Base.Threads
   include("threaded.jl")
end

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
       @sg_str,
       NucleotideSequence,
       RefSet, RefTx, RefGene,
       Multimap,
       AlignParam,
       ungapped_align,
       make_fqparser,
       isgzipped,
       load_refflat,
       load_gtf,
       fasta_to_index,
       process_reads!,
       process_paired_reads!,
       calculate_tpm!,
       rec_gene_em!,
       output_tpm,
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
