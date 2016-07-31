if VERSION >= v"0.5-"
   using Base.Test
else        
   using BaseTestNext 
   const Test = BaseTestNext            
end

using DataStructures
using BufferedStreams
using Bio.Seq
using FMIndexes
using IntArrays
using IntervalTrees
using Libz
using Distributions

include("../src/types.jl")
include("../src/timer.jl")
include("../src/bio_nuc_safepatch.jl")
include("../src/refset.jl")
include("../src/graph.jl")
include("../src/edges.jl")
include("../src/index.jl")
include("../src/align.jl")
include("../src/quant.jl")
include("../src/reads.jl")
include("../src/paired.jl")
include("../src/events.jl")
include("../src/io.jl")
include("../src/diff.jl")

@testset "Bio.Seq Patch" begin
   @test typeof(sg"GATGCA") == NucleotideSequence{SGNucleotide}
   @test reverse_complement(sg"GATGCA") == sg"TGCATC"
   @test reverse_complement(sg"LRS")    == sg"SRL" 
   
end
@testset "Splice Graphs" begin
   gtf = IOBuffer("# gtf file test
chr0\tTEST\texon\t6\t20\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t31\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t54\t62\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t76\t85\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t6\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t51\t62\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t76\t85\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t6\t20\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t31\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t54\t65\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t76\t90\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t11\t20\t.\t+\t.\tgene_id \"single\"; transcript_id \"ex1_single\";
chr0\tTEST\texon\t11\t20\t.\t-\t.\tgene_id \"single_rev\"; transcript_id \"single_rev\";
chr0\tTEST\texon\t11\t20\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"def_kiss\";
chr0\tTEST\texon\t21\t30\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"def_kiss\";
")

   flat = IOBuffer("# refflat file test (gtfToGenePred -genePredExt test.gtf test.flat)
def\tchr0\t+\t5\t85\t85\t85\t4\t5,30,53,75,\t20,40,62,85,\t0\tone\tnone\tnone\t-1,-1,-1,-1,
int1_alt3\tchr0\t+\t5\t85\t85\t85\t3\t5,50,75,\t40,62,85,\t0\tone\tnone\tnone\t-1,-1,-1,
apa_alt5\tchr0\t+\t5\t90\t90\t90\t4\t5,30,53,75,\t20,40,65,90,\t0\tone\tnone\tnone\t-1,-1,-1,-1,
ex1_single\tchr0\t+\t10\t20\t10\t20\t1\t10,\t20,\t0\tsingle\tnone\tnone\t-1,
")

   gtfref  = load_gtf( gtf )
   flatref = load_refflat( flat )

   @testset "Gene Annotation" begin

      for gene in keys(flatref)
         @test gtfref[gene].don    == flatref[gene].don
         @test gtfref[gene].acc    == flatref[gene].acc
         @test gtfref[gene].txst   == flatref[gene].txst
         @test gtfref[gene].txen   == flatref[gene].txen
         @test gtfref[gene].length == flatref[gene].length
      end

   end

                              # fwd     rev
   buffer1   = sg"AAAAA"      # 1-5     96-100
   utr5      = sg"TTATT"      # 6-10    91-95
   exon1     = sg"GCGGATTACA" # 11-20   81-90
   int1      = sg"TTTTTTTTTT" # 21-30   71-80
   exon2     = sg"GCATTAGAAG" # 31-40   61-70
   int2      = sg"GGGGGGGGGG" # 41-50   51-60
   exon3alt3 = sg"CCT"        # 51-53   48-50
   exon3def  = sg"CTATGCTAG"  # 54-62   39-47
   exon3alt5 = sg"TTC"        # 63-65   36-38
   int3      = sg"CCCCCCCCCC" # 66-75   26-35
   exon4     = sg"TTAGACAAGA" # 76-85   16-25
   apa       = sg"AATAA"      # 86-90   11-15
   buffer2   = sg"AAAAAAAAAA" # 91-100  1-10

   fwd = buffer1 * utr5 * exon1 * int1 * 
         exon2 * int2 * 
         exon3alt3 * exon3def * exon3alt5 * int3 * 
         exon4 * apa * buffer2
   rev = reverse_complement(fwd)

   genome = fwd * rev

   graphseq_one = sg"SL" * utr5 * exon1 * sg"LL" * 
               int1 * sg"RR" *
               exon2 * sg"LR" * 
               exon3alt3 * sg"RR" * 
               exon3def * sg"LL" * 
               exon3alt5 * sg"LR" *
               exon4 * sg"RS" * 
               apa * sg"RS"

   graphseq_sin = sg"SLGCGGATTACARS"
   graphseq_kis = sg"SLAAAAAAAAAALRTGTAATCCGCRS"

   graph_one = SpliceGraph( gtfref["one"], genome )
   graph_sin = SpliceGraph( gtfref["single"], genome )
   graph_rev = SpliceGraph( gtfref["single_rev"], genome)
   graph_kis = SpliceGraph( gtfref["kissing"], genome )

   println(graph_one)
   println(graph_sin)
   println(graph_rev)
   println(graph_kis)

   @testset "Graph Building" begin
      @test graph_one.seq == graphseq_one
      @test graph_sin.seq == graphseq_sin
      @test graph_kis.seq == graphseq_kis
   end

   @testset "Kmer Edges" begin

   end

   @testset "Alignment" begin

   end
   
end
