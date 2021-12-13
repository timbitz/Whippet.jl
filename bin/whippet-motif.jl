#!/usr/bin/env julia
# Tim Sterne-Weiler 2021

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

start = time_ns()
println( stderr, "Whippet $ver loading... " )

Pkg.activate(dir * "/..")

using Whippet
using Random
using ArgParse
using Glob
using BufferedStreams
using Libz
using BioSequences
using FASTX
using DataFrames

function parse_cmd()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--files"
      help     = "Pattern to glob.psi[.gz] (common-filename-segment [*.psi*]), or comma delimited list of filenames. ie. (-a sample_a) would work for sample_a-rep1.psi.gz,sample_a-rep2.psi.gz,..."
      arg_type = String
      default  = "*psi.gz"
    "--fasta"
      help     = "Genome sequence fasta file used with `--fasta` parameter when building the whippet index"
      arg_type = String
      required = true
    "--rbp-span"
      help     = "Size of cis-regulatory region upstream and downstream to score with RBNS-data (default:50)"
      arg_type = Int
      default  = 50

  end
  return parse_args(s)
end

function remove_ambiguous!( seq, replace_with::BioSequences.DNA=DNA_A )
   for (pos, nuc) in each(isambiguous, seq)
      seq[pos] = replace_with
   end
end

function read_fasta( filename::String )
   genome = Dict{String,ReferenceSequence}()
   reader = FASTA.Reader( bufferedinput( filename ) )

   for entry in reader
      println(stderr, "Processing $(identifier(entry))")
      seq = sequence(entry)
      remove_ambiguous!(seq)
      genome[identifier(entry)] = ReferenceSequence(seq)
   end

   genome
end

function parse_coord( coord::S ) where S <: AbstractString
   chr,inter = split(coord, ':')
   if '-' in inter
      left,right = map(x->parse(Int,x), split(inter, '-'))
      return chr,left,right
   else
      middle = parse(Int, inter)
      return chr,middle
   end
end

function cat_string( arr::Vector{A} ) where A
   str = "["
   for i in 1:length(arr)-1
      str *= string(arr[i])
      str *= ","
   end
   str *= string(arr[end])
   str *= "]"
   str
end

function score_motifs( ifile::String, 
                       ofile::String,
                       five_mat::Vector{Float64},
                       three_mat::DataFrame,
                       rbps::Vector{RNABindNSeq},
                       genome::Dict{String,ReferenceSequence},
                       rbpspan::Int,
                       gpsi::Bool=false )

   header    = true
   skipfive  = false
   skipthree = false

   used = Dict{String,Bool}()

   bufout = open(ofile, "w") |> ZlibDeflateOutputStream

   for l in eachline( bufferedinput(ifile) )

      if header
         rbpstr = join(map(x->x.name * "_" * x.conc, rbps), ',')
         tab_write(bufout, "Gene\tNode\tCoord\tStrand\tType")
         tab_write(bufout, "SpliceSite3p_UpstreamRBP:[$rbpstr]")
         tab_write(bufout, "SpliceSite3p_Seq")
         tab_write(bufout, "SpliceSite3p_MaxEnt")
         tab_write(bufout, "SpliceSite3p_DownstreamRBP:[$rbpstr]")
         tab_write(bufout, "SpliceSite5p_UpstreamRBP:[$rbpstr]")
         tab_write(bufout, "SpliceSite5p_Seq")
         tab_write(bufout, "SpliceSite5p_MaxEnt")
         end_write(bufout, "SpliceSite5p_DownstreamRBP:[$rbpstr]")
         header = false
         continue
      end

      spl  = split(l, '\t')
      pre,fval = gpsi ? (2,4) : (1,3)
      ci, si, ti = fval,fval+1,fval+2 #prefix, coord, strand, type

      key = join(spl[pre:ti], ':')
      if haskey(used, key)
         continue
      else
         used[key] = true
      end

      for i in pre:ti
         tab_write(bufout, spl[i])
      end

      if spl[ti] in ("XDI", "XDE", "AD")

         if spl[ti] == "AD"
            name,left,right = parse_coord(spl[ci])
            dcoord = spl[si][1] == '+' ? right : left
         else
            name,dcoord = parse_coord(spl[ci])
         end

         skipthree = true

      elseif spl[ti] in ("XAI", "XAE", "AA")

         if spl[ti] == "AA"
            name,left,right = parse_coord(spl[ci])
            acoord = spl[si][1] == '+' ? left : right
         else
            name,acoord = parse_coord(spl[ci])
         end

         skipfive = true


      elseif spl[ti] in ("RI", "XRI", "XEI")
         name,left,right = parse_coord(spl[ci])

         if spl[ti] == "RI" # Bug-adjust
            if spl[si][1] == '+'
               left = genome[name][left:left+1] == dna"GT" ? left : left + 1
            else
               left = genome[name][left:left+1] == dna"CT" ? left : left + 1
            end
         end

         dcoord = spl[si][1] == '+' ? left - 1 : right + 1
         acoord = spl[si][1] == '+' ? right + 1 : left - 1

      else
         name,left,right = parse_coord(spl[ci])

         dcoord = spl[si][1] == '+' ? right : left
         acoord = spl[si][1] == '+' ? left : right
      end

      if !skipthree
         if spl[si][1] == '+'
            thress = LongDNASeq(genome[name][acoord-20:acoord+2])
            upseq  = LongDNASeq(genome[name][acoord-rbpspan-6:acoord-6])
            doseq  = LongDNASeq(genome[name][acoord:acoord+rbpspan])
         else
            thress = LongDNASeq(genome[name][acoord-2:acoord+20])
            upseq  = LongDNASeq(genome[name][acoord-rbpspan:acoord])
            doseq  = LongDNASeq(genome[name][acoord+6:acoord+rbpspan])
            reverse_complement!(thress)
            reverse_complement!(upseq)
            reverse_complement!(doseq)
         end
         maxent = normalize_score( score_three(thress, three_mat) )
         uprbps = score_cis(upseq, rbps)
         dorbps = score_cis(doseq, rbps)

         tab_write(bufout, cat_string(uprbps))
         tab_write(bufout, string(thress))
         tab_write(bufout, string(maxent))
         tab_write(bufout, cat_string(dorbps))
      else
         tab_write(bufout, "NA\tNA\tNA\tNA")
      end

      if !skipfive
         if spl[si][1] == '+'
            fivess = LongDNASeq(genome[name][dcoord-2:dcoord+6])
            upseq  = LongDNASeq(genome[name][dcoord-rbpspan:dcoord])
            doseq  = LongDNASeq(genome[name][dcoord+6:dcoord+6+rbpspan])

         else
            fivess = LongDNASeq(genome[name][dcoord-6:dcoord+2])
            upseq  = LongDNASeq(genome[name][dcoord:dcoord+rbpspan])
            doseq  = LongDNASeq(genome[name][dcoord-rbpspan-6:dcoord-6])
            reverse_complement!(fivess)
            reverse_complement!(upseq)
            reverse_complement!(doseq)
         end
         maxent = normalize_score( score_five(fivess, five_mat) )
         uprbps = score_cis(upseq, rbps)
         dorbps = score_cis(doseq, rbps)

         tab_write(bufout, cat_string(uprbps))
         tab_write(bufout, string(fivess))
         tab_write(bufout, string(maxent))
         write(bufout, cat_string(dorbps))
      else
         write(bufout, "NA\tNA\tNA\tNA")
      end

      skipfive  = false
      skipthree = false

      write(bufout, "\n")
   end

   close(bufout)
end

function main()

   args  = parse_cmd()
   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )

   five_mat  = load_matrix5( dir * "/../motif/maxent/score5_matrix.txt.gz" )
   three_mat = load_matrix3( dir * "/../motif/maxent/score3_matrix.txt.gz" )
   rbps = map(load_rbns, glob("*.tsv.gz", fixpath(dir * "/../motif/5mers/")))

   genome   = read_fasta( args["fasta"] )
   rbpspan  = args["rbp-span"]

   name = basename(args["files"])

   ifiles = glob(name, replace(args["files"], "$name" => "")) 
   sfiles = map(x->(isgzipped(x) ? splitext(x)[1] : x), ifiles)
   ofiles = map(x->splitext(basename(x))[1] * ".motif.gz", sfiles)
   fexten = map(x->(split(x, '.')[end] == "gz" ? split(x, '.')[end-1] : split(x, '.')[end]), sfiles)

   Threads.@threads for i in 1:length(ifiles)
      if isfile(ofiles[i])
         continue
      end
      println(stderr, "Fetching motifs for $(ifiles[i])")
      score_motifs( ifiles[i], ofiles[i], five_mat, three_mat, rbps, genome, rbpspan, fexten[i] == "gpsi" )
   end

   println(stderr, "Whippet $ver done." )
end

@timer main()
