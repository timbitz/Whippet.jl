#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

start = time_ns()
println( stderr, "Whippet $ver loading... " )

Pkg.activate(dir * "/..")

using Whippet
using BioAlignments
using XAM
using Libz
using Serialization
using Nullables
using ArgParse

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--kmer", "-k"
      help = "Kmer size to use for exon-exon junctions (default 9)"
      arg_type = Int64
      default  = 9
    "--fasta"
      help = "File containg the genome in fasta, one entry per chromosome [.gz]"
      arg_type = String
      required = true
    "--gtf"
      help = "Gene anotation file in GTF format"
      arg_type = String
    "--bam"
      help = "Sorted and Indexed BAM file to supplement index with de novo splice-sites/exons/retained-introns. (Must have corresponding .bai file)"
      arg_type = String
    "--index", "-x"
      help = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = String
      default  = "$dir/../index/graph"
    "--suppress-low-tsl"
      help = "Ignore low quality transcript annotations with TSL2+"
      action   = :store_true
    "--bam-min-reads"
      help = "Minimum number of reads supporting a splice-site from BAM to include in index."
      arg_type = Int
      default  = 1
    "--bam-both-novel"
      help = "Allow spliced reads from BAM where neither of the splice-sites are known (default off)."
      action   = :store_true
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()

   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds." )

   if args["gtf"] == nothing
      println(stderr, "ERROR: Must supply gene annotation file using `--gtf`!!")
      exit(1)
   else
      annotype = "gtf"
      annotxt  = "GTF"
   end

   flat = fixpath( args[annotype] )
   println(stderr, "Loading $annotxt file: $flat")
   fh = open( flat , "r")
   if isgzipped( flat )
      fh = fh |> x->ZlibInflateInputStream(x, reset_on_end=true)
   end

   if args["bam"] != nothing
      bam = fixpath( args["bam"] )
      isfile(bam) || error("ERROR: --bam parameter used, but cannot find .bam file at $bam !")
      isfile(bam * ".bai") || error("ERROR: --bam parameter used, but no .bai index found for .bam file! Cannot set-up random access to $bam !")
      println(stderr, "Loading BAM file for random-access: $bam")
      bamreadr = open(BAM.Reader, bam, index=bam * ".bai")
      @timer ref = load_gtf(fh, suppress=args["suppress-low-tsl"],
                                usebam=true,
                                bamreader=Nullable(bamreadr),
                                bamreads=args["bam-min-reads"],
                                bamoneknown=!args["bam-both-novel"])
   else
      @timer ref = load_gtf(fh, suppress=args["suppress-low-tsl"])
   end

   println(stderr, "Indexing transcriptome...")
   @timer graphome = fasta_to_index( fixpath( args["fasta"] ), ref, kmer=args["kmer"] )

   println(stderr, "Serializing splice graph index...")
   indexname = hasextension(args["index"], ".jls") ? args["index"] : args["index"] * ".jls"
   open(indexname, "w") do io
      @timer serialize( io, graphome )
   end

   println(stderr, "Printing a map of Whippet nodes to putative full exons... $(args["index"] * ".exons.tab.gz")")
   output_exons( args["index"] * ".exons.tab.gz", graphome )

end

@timer main()
