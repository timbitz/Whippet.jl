#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

tic()
println( STDERR, "Whippet $ver loading and compiling... " )

using ArgParse

push!( LOAD_PATH, dir * "/../src" )
using Whippet
using Libz

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
    "--flat"
      help = "Gene annotation file in RefFlat format"
      arg_type = String
    "--gtf"
      help = "Gene anotation file in GTF format"
      arg_type = String
    "--index", "-x"
      help = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph)"
      arg_type = String
      default = "$dir/../index/graph"
    "--suppress-low-tsl"
      help = "Ignore low quality transcript annotations with TSL2+"
      action   = :store_true
  end
  return parse_args(s)
end

function main()

   args = parse_cmd()
   
   println(STDERR, " $( round( toq(), 6 ) ) seconds." )

   if args["gtf"] == nothing && args["flat"] == nothing
      println(STDERR, "ERROR: Must supply gene annotation file using `--gtf` or `--flat`!!")
      exit(1)
   elseif args["gtf"] != nothing
      annotype = "gtf"
      annotxt  = "GTF"
   else
      annotype = "flat"
      annotxt  = "Refflat"
   end

   println(STDERR, "Loading $annotxt file...")
   flat = fixpath( args[annotype] )
   fh = open( flat , "r")
   if isgzipped( flat )
      fh = fh |> x->ZlibInflateInputStream(x, reset_on_end=true)
   end
   @timer ref = annotype == "gtf" ? load_gtf(fh, suppress=args["suppress-low-tsl"]) : load_refflat(fh)

   println(STDERR, "Indexing transcriptome...")
   @timer graphome = fasta_to_index( fixpath( args["fasta"] ), ref, kmer=args["kmer"] )

   #=
   println(STDERR, "Saving Annotations...")
   open("$(args["index"])_anno.jls", "w") do fh
      @timer serialize(fh, ref)
   end
   =#

   println(STDERR, "Serializing splice graph index...")
   indexname = hasextension(args["index"], ".jls") ? args["index"] : args["index"] * ".jls"
   open(indexname, "w") do io
      @timer serialize( io, graphome )
   end

end

@timer main()
