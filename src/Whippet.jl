__precompile__()

module Whippet

import BioSymbols
import BioSequences
import BioAlignments
import FASTX
import XAM

import DataStructures: SortedSet

using StatsBase: midpoints
using BioSymbols
using BioSequences
using BioAlignments
using GenomicFeatures
using FASTX
using XAM
using BufferedStreams
using Distributions
using FMIndexes
using IntArrays
using IntervalTrees
using Libz
using LinearAlgebra
using Nullables
using Printf
using Random
using CSV
using DataFrames

include("types.jl")
include("timer.jl")
include("sgkmer.jl")
include("fmindex_patch.jl")
include("bam.jl")
include("refset.jl")
include("graph.jl")
include("bias.jl")
include("edges.jl")
include("index.jl")
include("record.jl")
include("align.jl")
include("quant.jl")
include("reads.jl")
include("paired.jl")
include("events.jl")
include("io.jl")
include("diff.jl")
include("motif.jl")

export SpliceGraph,
       SpliceGraphQuant,
       GraphLib,
       GraphLibQuant,
       GeneInfo,
       TxInfo,
       SGNode,
       SGNodeMeta,
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
       SGAlignSingle,
       SGAlignPaired,
       @sg_str,
       reverse_complement,
       kmer,
       NucleotideSequence,
       RefSet, RefTx, RefGene,
       MultiMapping,
       DefaultBiasMod,
       DefaultCounter,
       PrimerBiasMod,
       PrimerBiasCounter,
       GCBiasMod,
       GCBiasCounter,
       JointBiasMod,
       JointBiasCounter,
       AlignParam,
       ungapped_align,
       make_fqparser,
       make_collated_parser,
       make_single_parser,
       read_collated_bam,
       check_samtools,
       ident_to_fastq_url,
       hasextension,
       isgzipped,
       ispaired_bamfile,
       load_refflat,
       load_gtf,
       fasta_to_index,
       checkversion,
       process_reads!,
       process_paired_reads!,
       build_equivalence_classes!,
       set_global_readcount,
       total_multi,
       primer_adjust!,
       gc_adjust!,
       primer_normalize!,
       gc_normalize!,
       calculate_tpm!,
       gene_em!,
       set_gene_tpm!,
       output_tpm,
       output_stats,
       output_junctions,
       output_exons,
       assign_ambig!,
       effective_lengths!,
       global_bias,
       process_events,
       fixpath,
       @timer,
       timer_print,
       open_streams,
       process_psi_files,
       tab_write,
       RNABindNSeq,
       load_rbns,
       load_matrix5,
       load_matrix3,
       score_cis,
       score_five,
       score_three,
       normalize_score
end
