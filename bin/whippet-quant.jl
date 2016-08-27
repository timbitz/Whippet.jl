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
      default  = fixpath( "./output" )
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
      help     = "Allowable number of mismatches in alignment (counted as 1-10^(-phred/10))"
      arg_type = Int
      default  = 3
    "--score-min", "-S"
      help     = "Minimum percent matching (matches - mismatches) / read_length"
      arg_type = Float64
      default  = 0.6
    "--psi-body-read"
      help     = "Allow exon-body reads in quantification of PSI values"
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
    "--curl"
      help     = "FASTQ files are URLs to download/process on the fly"
      action   = :store_true
    "--ebi"
      help     = "Retrieve FASTQ files from ebi.ac.uk using seq run id (ie. SRR1199003). (sets --curl=true)"
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

   # do we need to fetch data from ebi.ac.uk?
   if args["ebi"]
      ebi_res = ident_to_fastq_url( args["filename.fastq[.gz]"] )
      ispaired = ebi_res.paired
      args["curl"] = true
      args["filename.fastq[.gz]"] = ebi_res.fastq_1_url
      if ispaired
         args["paired_mate.fastq[.gz]"] = ebi_res.fastq_2_url
      end
      ebi_res.success || error("Could not fetch data from ebi.ac.uk!!")
   end

   const parser,iobuf,rref = args["curl"] ? make_http_fqparser( args["filename.fastq[.gz]"],
                                                           encoding=enc, forcegzip=args["force-gz"] ) : 
                                            make_fqparser( fixpath(args["filename.fastq[.gz]"]), 
                                                           encoding=enc, forcegzip=args["force-gz"] )
   if ispaired
      const mate_parser,mate_iobuf,mate_rref = args["curl"] ? make_http_fqparser( args["paired_mate.fastq[.gz]"],
                                                                             encoding=enc, forcegzip=args["force-gz"] ) :
                                                              make_fqparser( fixpath(args["paired_mate.fastq[.gz]"]), 
                                                                             encoding=enc, forcegzip=args["force-gz"] )
   end

   if nprocs() > 1
      println(STDERR, "Whippet does not currrently support nprocs() > 1")
      return #TODO: first implementation was too slow, ie too much communication overhead
   else
      println(STDERR, "Processing reads...")
      if ispaired
         @timer mapped,total,readlen = process_paired_reads!( parser, mate_parser, param, lib, quant, multi, 
                                                              sam=args["sam"], qualoffset=enc_offset, 
                                                              iobuf=iobuf, mate_iobuf=mate_iobuf,
                                                              rref=rref, mate_rref=mate_rref, http=args["curl"] )
         readlen = round(Int, readlen)
         println(STDERR, "Finished mapping $mapped paired-end reads of length $readlen each out of a total $total mate-pairs...")
      else
         @timer mapped,total,readlen = process_reads!( parser, param, lib, quant, multi, 
                                                       sam=args["sam"], qualoffset=enc_offset,
                                                       iobuf=iobuf, rref=rref, http=args["curl"] )
         readlen = round(Int, readlen)
         println(STDERR, "Finished mapping $mapped single-end reads of length $readlen out of a total $total reads...")
      end
   end

   # TPM_EM
   println(STDERR, "Calculating expression values and MLE for $( length(multi) ) repetitive reads with EM...")
   calculate_tpm!( quant, readlen=readlen )
   @timer iter = rec_gene_em!( quant, multi, sig=6, readlen=readlen, max=1000 ) 
   println(STDERR, "Finished calculating transcripts per million (TpM) after $iter iterations of EM...")

   output_tpm( args["out"] * ".tpm.gz", lib, quant )
   output_stats( args["out"] * ".map.gz", lib, quant, param, total, mapped, length(multi), readlen )

   println(STDERR, "Assigning multi-mapping reads based on maximum likelihood estimate..")
   # Now assign multi to edges.
   @timer assign_ambig!( quant, multi, ispaired=ispaired )

   println(STDERR, "Calculating effective lengths...")
   @timer effective_lengths!( lib, quant, readlen - 19, 9-1) #min(readlen - param.score_min, 9-1) )
   @timer bias_ave,bias_var = global_bias( quant )
   println(STDERR, "Global bias is $bias_ave +/- $bias_var ")
   println(STDERR, "Calculating maximum likelihood estimate of events..." )
   @timer process_events( args["out"] * ".psi.gz" , lib, quant, isnodeok=args["psi-body-read"], iscircok=args["circ"] )
   println(STDERR, "Whippet $ver done." )
end

main()
