#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

tic()
println( STDERR, "Whippet $ver loading and compiling... " )

using ArgParse

push!( LOAD_PATH, dir * "/../src" )
import SpliceGraphs
using SpliceGraphs

function parse_cmd()
  s = ArgParseSettings()
  # TODO finish options...
  @add_arg_table s begin
    "filename.fastq[.gz]"
      arg_type = ASCIIString
      required = true
    "paired_mate.fastq[.gz]"
      arg_type = ASCIIString
      required = false
    "--index", "-x"
      help     = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = ASCIIString
      default  = fixpath( "$(dir)/../index/graph" )
    "--out", "-o"
      help     = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = ASCIIString
      default  = fixpath( "$(dir)/../output" )
    "--sam", "-s"
      help     = "Should SAM format be sent to stdout?"
      action   = :store_true
    "--seed-len", "-L"
      help     = "Seed length"
      arg_type = Int
      default  = 18
    "--seed-try", "-M"
      help     = "Number of failed seeds to try before giving up"
      arg_type = Int
      default  = 3
    "--seed-tol", "-T"
      help     = "Number of seed hits to tolerate"
      arg_type = Int
      default  = 4
    "--seed-buf", "-B"
      help     = "Ignore this many bases from beginning and end of read for seed"
      arg_type = Int
      default  = 5
    "--seed-inc", "-I"
      help     = "Number of bases to increment seed each iteration"
      arg_type = Int
      default  = 18
    "--pair-range", "-P"
      help     = "Seeds for paired end reads must match within _ bases of one another"
      arg_type = Int
      default  = 2500
    "--mismatches", "-X"
      help     = "Allowable number of mismatches in alignment"
      arg_type = Int
      default  = 3
    "--score-min", "-S"
      help     = "Minimum alignment score (matches - mismatches)"
      arg_type = Int
      default  = 45
    "--junc-only", "-j"
      help     = "Only use junction reads, no internal exon reads will be considered."
      action   = :store_true
    "--stranded"
      help     = "Is the data strand specific? If so, increase speed with this flag"
      action   = :store_true
    "--rev-pair"
      help     = "Is the second mate the reverse complement of the first? If so, increase speed with this flag"
      action   = :store_true
    "--phred-33" 
      help     = "Qual string is encoded in Phred+33 integers (default)"
      action   = :store_true
    "--phred-64"
      help     = "Qual string is encoded in Phred+64 integers"
      action   = :store_true
    "--no-circ"
      help     = "Do not allow back/circular splicing"
      action   = :store_false
    "--no-tpm"
      help     = "Should tpm file be sent to output/prefix.tpm.gz? (default on)"
      action   = :store_true
    "--force-gz"
      help     = "Regardless of suffix, consider read input as gzipped"
      action   = :store_true
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(STDERR, " $( round( toq(), 6 ) ) seconds" )

   indexpath = fixpath( args["index"] )
   println(STDERR, "Loading splice graph index... $( indexpath ).jls")
   @timer const lib = open(deserialize, "$( indexpath ).jls")

   println(STDERR, "Loading annotation index... $( indexpath )_anno.jls")
   @timer const anno = open(deserialize, "$( indexpath )_anno.jls")

   const ispaired = args["paired_mate.fastq[.gz]"] != nothing ? true : false

   const param = AlignParam( args, ispaired, kmer=lib.kmer ) 
   const quant = GraphLibQuant( lib, anno )
   const multi = Vector{Multimap}()

   const enc        = args["phred-64"] ? Bio.Seq.ILLUMINA15_QUAL_ENCODING : Bio.Seq.ILLUMINA18_QUAL_ENCODING
   const enc_offset = args["phred-64"] ? 64 : 33

   const parser = make_fqparser( fixpath(args["filename.fastq[.gz]"]), encoding=enc, forcegzip=args["force-gz"] )
   if ispaired
      const mate_parser = make_fqparser( fixpath(args["paired_mate.fastq[.gz]"]), encoding=enc, forcegzip=args["force-gz"] )
   end

   if nprocs() > 1
      #include("align_parallel.jl")
      # Load Fastq files in chunks
      # Parallel reduction loop through fastq chunks
      println(STDERR, "Whippet does not currrently support nprocs() > 1")
      return #TODO: first implementation was too slow, ie too much communication overhead
   else
      println(STDERR, "Processing reads...")
      if ispaired
         @timer mapped,total,readlen = process_paired_reads!( parser, mate_parser, param, lib, quant, multi, 
                                                              sam=args["sam"], qualoffset=enc_offset )
         readlen = round(Int, readlen)
         println(STDERR, "Finished mapping $mapped paired-end reads of length $readlen each out of a total $total mate-pairs...")
      else
         @timer mapped,total,readlen = process_reads!( parser, param, lib, quant, multi, 
                                                       sam=args["sam"], qualoffset=enc_offset )
         readlen = round(Int, readlen)
         println(STDERR, "Finished mapping $mapped single-end reads of length $readlen out of a total $total reads...")
      end
   end

   # TPM_EM
   println(STDERR, "Calculating expression values and MLE for $( length(multi) ) repetitive reads with EM...")
   calculate_tpm!( quant, readlen=readlen )
   @timer iter = rec_gene_em!( quant, multi, sig=6, readlen=readlen, max=1000 ) 
   println(STDERR, "Finished calculating transcripts per million (TpM) after $iter iterations of EM...")

   if !args["no-tpm"]
      output_tpm( args["out"] * ".tpm.gz", lib, quant )
   end

   println(STDERR, "Assigning multi-mapping reads based on maximum likelihood estimate..")
   # Now assign multi to edges.
   @timer assign_ambig!( quant, multi, ispaired=ispaired )

   println(STDERR, "Calculating effective lengths...")
   @timer effective_lengths!( lib, quant, readlen - 19, min(readlen - param.score_min, 9-1) )
   @timer bias_ave,bias_var = global_bias( quant )
   println(STDERR, "Global bias is $bias_ave +/- $bias_var ")
   println(STDERR, "Calculating maximum likelihood estimate of events..." )
   @timer process_events( args["out"] * ".psi.gz" , lib, quant, isnodeok=!args["junc-only"] )
   println(STDERR, "Whippet $ver done." )
end

main()
