#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = readline(open(dir * "/VERSION"))

start = time_ns()
println( stderr, "Whippet $ver loading... " )

Pkg.activate(dir * "/..")

using Whippet
using Serialization
using ArgParse

function parse_cmd()
  s = ArgParseSettings()
  # TODO finish options...
  @add_arg_table s begin
    "filename.fastq[.gz]"
      arg_type = String
      required = true
    "paired_mate.fastq[.gz]"
      arg_type = String
      required = false
    "--index", "-x"
      help     = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph[.jls])"
      arg_type = String
      default  = fixpath( "$dir/../index/graph" )
    "--out", "-o"
      help     = "Where should the gzipped output go 'dir/prefix'?"
      arg_type = String
      default  = fixpath( "./output" )
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
      default  = 5000
    "--mismatches", "-X"
      help     = "Allowable number of mismatches in alignment (counted as 1-10^(-phred/10))"
      arg_type = Int
      default  = 3
    "--score-min", "-S"
      help     = "Minimum percent matching (matches - mismatches) / read_length"
      arg_type = Float64
      default  = 0.7
    "--sam"
      help     = "Should SAM format be sent to stdout?"
      action   = :store_true
    "--biascorrect"
      help     = "Apply fragment GC-content and 5' Sequence bias correction methods for more stable PSI values"
      action   = :store_true
    "--stranded"
      help     = "Is the data strand specific in fwd orientation? If so, increase speed with this flag"
      action   = :store_true
    "--pair-same-strand"
      help     = "Whippet by default tries to align fwd/rev pairs, if your data is fwd/fwd or rev/rev set this flag"
      action   = :store_true
    "--phred-33"
      help     = "Qual string is encoded in Phred+33 integers (default)"
      action   = :store_true
    "--phred-64"
      help     = "Qual string is encoded in Phred+64 integers"
      action   = :store_true
    "--circ"
      help     = "Allow back/circular splicing, this will allow output of `BS`-type lines"
      action   = :store_true
    "--force-gz"
      help     = "Regardless of suffix, consider read input as gzipped"
      action   = :store_true
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds." )

   indexpath = fixpath( args["index"] )
   indexname = hasextension(indexpath, ".jls") ? indexpath : indexpath * ".jls"
   println(stderr, "Loading splice graph index... $indexname")
   @timer lib = open(deserialize, indexname)

   ispaired = args["paired_mate.fastq[.gz]"] != nothing

    args["filename.fastq[.gz]"] = fixpath(args["filename.fastq[.gz]"])
    if ispaired
       args["paired_mate.fastq[.gz]"] = fixpath(args["paired_mate.fastq[.gz]"])
    end

   ContainerType = ispaired ? SGAlignPaired : SGAlignSingle

   if args["biascorrect"]
      CounterType   = JointBiasCounter
      BiasModelType = JointBiasMod
   else
      CounterType   = DefaultCounter
      BiasModelType = DefaultBiasMod
   end

   param = AlignParam( args, ispaired, kmer=lib.kmer )
   quant = GraphLibQuant{ContainerType,CounterType}( lib )
   multi = MultiMapping{ContainerType,CounterType}()
   mod   = BiasModelType()

   enc_offset = args["phred-64"] ? 64 : 33

   parser = make_fqparser( args["filename.fastq[.gz]"], forcegzip=args["force-gz"] )

   if ispaired
      mate_parser = make_fqparser( args["paired_mate.fastq[.gz]"], forcegzip=args["force-gz"] )
   end

   #TODO: first implementation was too slow, ie too much communication overhead
   # try speeding up by using Distributed.nprocs > 1
   println(stderr, "Processing reads from file...")
    if ispaired
       println(stderr, "FASTQ_1: " * args["filename.fastq[.gz]"])
       println(stderr, "FASTQ_2: " * args["paired_mate.fastq[.gz]"])
       @timer mapped,totreads,readlen = process_paired_reads!( parser, mate_parser, param, lib, quant, multi, mod,
                                                           sam=args["sam"], qualoffset=enc_offset )
       readlen = Int(floor(readlen))
       println(stderr, "Finished mapping $mapped paired-end reads of length $readlen each out of a total of $totreads mate-pairs...")
    else
       println(stderr, "FASTQ: " * args["filename.fastq[.gz]"])
       @timer mapped,totreads,readlen = process_reads!( parser, param, lib, quant, multi, mod,
                                                     sam=args["sam"], qualoffset=enc_offset )
       readlen = Int(floor(readlen))
       println(stderr, "Finished mapping $mapped single-end reads of length $readlen out of a total of $totreads reads...")
    end

   # TPM_EM
   println(stderr, "Calculating expression values and MLE of equivalence classes with EM:")
   primer_normalize!( mod )
   primer_adjust!( quant, mod )
   primer_adjust!( multi, mod )
   println(stderr, "- $( length(multi.map) ) multi-gene mapping read equivalence classes...")
   build_equivalence_classes!( quant, lib, assign_long=true )
   println(stderr, "- $( length(quant.classes) ) multi-isoform equivalence classes...")
   calculate_tpm!( quant, readlen=readlen )
   @timer iter = gene_em!( quant, multi, sig=1, readlen=readlen, maxit=10000 )
   println(stderr, "Finished calculating transcripts per million (TpM) after $iter iterations of EM...")
   set_gene_tpm!( quant, lib )
   gc_normalize!( mod, lib, quant )
   gc_adjust!( quant, mod )
   gc_adjust!( multi, mod )

   output_tpm( args["out"], lib, quant )
   output_stats( args["out"] * ".map.gz", lib, quant, param, indexpath, totreads, mapped, Int(round(total_multi(multi))), readlen, ver )
   output_junctions( args["out"] * ".jnc.gz", lib, quant )

   println(stderr, "Assigning multi-mapping reads based on maximum likelihood estimate..")
   # Now assign multi to edges.
   @timer assign_ambig!( quant, lib, multi )

   effective_lengths!( lib, quant, 1, 0)
   println(stderr, "Calculating maximum likelihood estimate of events..." )
   @timer process_events( args["out"] * ".psi.gz" , lib, quant, isnodeok=false, iscircok=args["circ"], readlen=readlen )
   println(stderr, "Whippet $ver done." )
end

@timer main()
