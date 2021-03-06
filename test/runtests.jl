using Test
using Serialization

import BioSymbols
import BioSequences
import BioAlignments
import XAM
import FASTX
using BioSymbols
using BioSequences
using BioAlignments
using GenomicFeatures
using XAM
using FASTX
using BufferedStreams
using DataStructures
using Distributions
using FMIndexes
using IntArrays
using InteractiveUtils
using IntervalTrees
using Libz
using LinearAlgebra
using Nullables
using Printf
using Random

include("../src/types.jl")
include("../src/timer.jl")
include("../src/sgkmer.jl")
include("../src/fmindex_patch.jl")
include("../src/bam.jl")
include("../src/refset.jl")
include("../src/graph.jl")
include("../src/bias.jl")
include("../src/edges.jl")
include("../src/index.jl")
include("../src/record.jl")
include("../src/align.jl")
include("../src/quant.jl")
include("../src/reads.jl")
include("../src/paired.jl")
include("../src/events.jl")
include("../src/io.jl")
include("../src/diff.jl")

@testset "SG Kmers & Seq" begin
   @test sgkmer(dna"ATG") == BioSequences.DNAMer(dna"ATG")
   @test sgkmer("ATG") == BioSequences.DNAMer(dna"ATG")
   @test isa(sgkmer(dna"ATG"), SGKmer{3})
   @test kmer_index( sgkmer(dna"ATG") ) == kmer_index( BioSequences.DNAMer(dna"ATG") )
   @test Int(kmer_index_trailing( dna"ATG" ))+1 == kmer_index( BioSequences.DNAMer(dna"ATG") )
   @test Int(kmer_index_trailing( dna"ATGR" )) == 0
   seq = dna""
   for i in 1:32
      seq = i % 2 == 0 ? seq * dna"C" : seq * dna"T"
      curkmer = sgkmer(seq)
      @test isa( curkmer, SGKmer{i} )
   end
end

@testset "Bias Models" begin
   # default 'empty' model
   def = zero(DefaultCounter)
   @test get(def) == 0.0
   def = one(DefaultCounter)
   @test get(def) == 1.0
   push!( def, 1.0 )
   @test get(def) == 2.0
   gc_adjust!( def, DefaultBiasMod() )
   primer_adjust!( def, DefaultBiasMod() )

   function randdna(n)
      return BioSequence{DNAAlphabet{2}}(rand([DNA_A, DNA_C, DNA_G, DNA_T], n))
   end

   function test_uniform_bias!( mod::B ) where B <: BiasModel
      for i in 1:5000000
         primer_count!( mod, randdna(36) )
      end
      primer_normalize!( mod )
      for i in 1:length(mod.fore)
         @test 0.85 < mod.back[i] / mod.fore[i] < 1.15
      end
   end

   # primer count model
   mod = JointBiasMod( 5 )
   test_uniform_bias!( mod )
   cnt = JointBiasCounter()
   @test cnt.count == 0.0
   hept = kmer_index_trailing(UInt16, dna"ATGAC")
   push!( cnt, hept )
   push!( cnt, hept )
   adjust!( cnt, mod )
   @test get(cnt) == value!( mod, hept )

   # gc content bias
   seq = dna"AAAAAGCGCG" ^ 5
   egc = ExpectedGC( seq )
   @test length(egc) == 20
   @test egc[ Int(div(gc_content(seq), 0.05)+1) ] == 1.0
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
chr0\tTEST\texon\t11\t30\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"ret_kiss\";
")

#=   flat = IOBuffer("# refflat file test (gtfToGenePred -genePredExt test.gtf test.flat)
def\tchr0\t+\t5\t85\t85\t85\t4\t5,30,53,75,\t20,40,62,85,\t0\tone\tnone\tnone\t-1,-1,-1,-1,
int1_alt3\tchr0\t+\t5\t85\t85\t85\t3\t5,50,75,\t40,62,85,\t0\tone\tnone\tnone\t-1,-1,-1,
apa_alt5\tchr0\t+\t5\t90\t90\t90\t4\t5,30,53,75,\t20,40,65,90,\t0\tone\tnone\tnone\t-1,-1,-1,-1,
ex1_single\tchr0\t+\t10\t20\t10\t20\t1\t10,\t20,\t0\tsingle\tnone\tnone\t-1,
")=#

   gtfref  = load_gtf( gtf )
   #flatref = load_refflat( flat )

   @testset "Gene Annotations" begin

      @test gtfref["one"].don == (20,40,62,65)
      @test gtfref["one"].acc == (31,51,54,76)
      @test gtfref["one"].noveldon == tuple()
      @test gtfref["one"].novelacc == tuple()

   end

   gtf.ptr = 1
   bamreadr = open(BAM.Reader, "test.sort.bam", index="test.sort.bam.bai" )
   supgtfref = load_gtf( gtf, usebam=true, bamreader=Nullable(bamreadr), bamreads=2 )

   @testset "BAM Supplemented Annotation" begin

      @test supgtfref["one"].don == (20,40,45,62,65)
      @test supgtfref["one"].acc == (31,51,54,76)
      @test supgtfref["one"].noveldon == tuple(45)
      @test supgtfref["one"].novelacc == tuple()

   end

                               # fwd     rev
   buffer1   = dna"AAAAA"      # 1-5     96-100
   utr5      = dna"TTATT"      # 6-10    91-95
   exon1     = dna"GCGGATTACA" # 11-20   81-90
   int1      = dna"TTTTTTTTTT" # 21-30   71-80
   exon2     = dna"GCATTAGAAG" # 31-40   61-70
   int2      = dna"GGGGGGGGGG" # 41-50   51-60
   exon3alt3 = dna"CCT"        # 51-53   48-50
   exon3def  = dna"CTATGCTAG"  # 54-62   39-47
   exon3alt5 = dna"TTC"        # 63-65   36-38
   int3      = dna"CCCCCCCCCC" # 66-75   26-35
   exon4     = dna"TTAGACAAGA" # 76-85   16-25
   apa       = dna"AATAA"      # 86-90   11-15
   buffer2   = dna"AAAAAAAAAA" # 91-100  1-10

   fwd = buffer1 * utr5 * exon1 * int1 *
         exon2 * int2 *
         exon3alt3 * exon3def * exon3alt5 * int3 *
         exon4 * apa * buffer2
   rev = reverse_complement(fwd)

   genome = fwd * rev

   expected_one = utr5 * exon1 *
               int1 *
               exon2 *
               exon3alt3 *
               exon3def *
               exon3alt5 *
               exon4 *
               apa

   expected_sin = dna"GCGGATTACA"
   expected_kis = dna"AAAAAAAAAATGTAATCCGC"

   kmer_size = 2 # good test size

   graph_one = SpliceGraph( gtfref["one"], genome, kmer_size )
   graph_sin = SpliceGraph( gtfref["single"], genome, kmer_size )
   graph_rev = SpliceGraph( gtfref["single_rev"], genome, kmer_size )
   graph_kis = SpliceGraph( gtfref["kissing"], genome, kmer_size )

   @testset "Graph Building" begin
      @test graph_one.seq == expected_one
      @test graph_sin.seq == expected_sin
      @test graph_kis.seq == expected_kis

      @test length(graph_one.annopath) == length(gtfref["one"].reftx)
      @test graph_one.annopath[1] == BitSet([1,3,5,7]) # path of def
      @test graph_one.annopath[2] == BitSet([1,2,3,4,5,7]) # path of int1_alt3
      @test graph_one.annopath[3] == BitSet([1,3,5,6,7,8]) # path of apa_alt5
   end

   # Build Index (from index.jl)
   xcript  = dna""
   xoffset = Vector{UInt64}()
   xgenes  = Vector{GeneName}()
   xinfo   = Vector{GeneInfo}()
   xlength = Vector{Float64}()
   xgraph  = Vector{SpliceGraph}()

   runoffset = 0

   for g in keys(gtfref)
      println(stderr, g)
      curgraph = SpliceGraph( gtfref[g], genome, kmer_size )
      xcript  *= curgraph.seq
      push!(xgraph, curgraph)
      push!(xgenes, g)
      push!(xinfo, gtfref[g].info )
      push!(xlength, gtfref[g].length )
      push!(xoffset, runoffset)
      runoffset += length(curgraph.seq)
   end

   fm = FMIndex(twobit_enc(xcript), 4, r=1, program=:SuffixArrays, mmap=true)

   edges = build_edges( xgraph, kmer_size )

   #println(edges.left)
   #println(edges.right)

   lib = GraphLib( xoffset, xgenes, xinfo, xlength, xgraph, edges, fm, true, kmer_size )

   @testset "Kmer Edges" begin
      left  = [dna"CA", dna"AG", dna"AG", dna"TC", dna"AA"]
      right = [dna"GC", dna"CC", dna"CT", dna"TT", dna"TG"]
      lkmer = map( x->kmer_index(SGKmer(x)), left )
      rkmer = map( x->kmer_index(SGKmer(x)), right )
      #println(lkmer[1])
      #println(rkmer[1])
      for i in 1:4^kmer_size
         if i in lkmer
            @test isassigned(edges.left, i)
            @test typeof(edges.left[i]) == Vector{SGNode}
            @test issorted(edges.left[i], lt=sortlt)
         else
            @test !isassigned(edges.left, i)
         end
         if i in rkmer
            @test isassigned(edges.right, i)
            @test typeof(edges.right[i]) == Vector{SGNode}
            @test issorted(edges.right[i], lt=sortlt)
         else
            @test !isassigned(edges.right, i)
         end
      end
      exon1_lind = lkmer[1]
      exon2_rind = rkmer[1]
      @test intersect( edges.left[exon1_lind], edges.right[exon2_rind] ) == edges.right[exon2_rind]
      @test edges.left[exon1_lind] ∩ edges.right[exon2_rind] == edges.right[exon2_rind]
   end

   @testset "Saving and Loading Index" begin
      println(stderr, "Saving test index...")
      println(stderr, lib)
      @timer open("test_index.jls", "w+") do io
         serialize( io, lib )
      end
      println(stderr, "Loading test index...")
      @timer lib = open(deserialize, "test_index.jls")
      println(stderr, lib)
   end

   # store general data structures outside of testset
   score_range = 0.05
   param = AlignParam( 1, 2, 4, 4, 4, 5, 1, 2, 1000, score_range, 0.7,
                          false, false, true, false, true )
   println(map( x->length(x.annoname), lib.graphs ))
   gquant = GraphLibQuant{SGAlignSingle,DefaultCounter}( lib )

   @testset "Alignment, SAM Format, Equivalence Classes" begin
      # reads
      fastq = IOBuffer("@1S10M%11,20%exon1
NGCGGATTACA
+
#BBBBBBBBBB
@1S9M%54,62%exon3def
NCTATGCTAG
+
#BBBBBBBBB
@1S15M%51,65%alt3-exon3-alt5
NCCTCTATGCTAGTTC
+
#BBBBBBBBBBBBBBB
@11M10N10M%10,40%exon1-exon2
NGCGGATTACAGCATTAGAAG
+
#BBBBBBBBBBBBBBBBBBBB
@5M10N6M%16,36%exon1trunc-exon2trunc
TTACAGCATTN
+
BBBBBBBBBB#
@10M33N9M1S%11,62%exon1-exon3def
GCGGATTACACTATGCTAGN
+
BBBBBBBBBBBBBBBBBBB#
@10M33N9M%11,62%exon1-exon3def:rc
CTAGCATAGTGTAATCCGC
+
BBBBBBBBBBBBBBBBBBB
@10M55N10M1S%11,85%exon1-exon4full
GCGGATTACATTAGACAAGAN
+
BBBBBBBBBBBBBBBBBBBB#
@11M55N4M%10,78%exon1-exon4_4bp
NGCGGATTACATTAG
+
#IIIIIIIIIIIIII
@2M55N10M%19,85%exon1_2bp-exon4:rc
TCTTGTCTAATG
+
IIIIIIIIIIII
")

      DNASeqType = BioSequences.BioSequence{BioSequences.DNAAlphabet{2}}
      fqparse = FASTQ.Reader( BufferedInputStream(fastq), fill_ambiguous=DNA_C )
      reads  = allocate_fastq_records( 10 )
      read_chunk!( reads, fqparse )

      @test length(reads) == 10
      for r in reads
         fill!( r, 33 )
         println(stderr, r)
         readname = string(r)

         align = ungapped_align( param, lib, r )
         println(stderr, align)

         flush(stderr)
         @test !isnull( align )
         @test length( align.value ) >= 1
         @test all(map( x->x.isvalid, align.value))
         scores = map( x->identity(x, length(r.sequence)), align.value )
         @test maximum(scores) - minimum(scores) <= score_range

         if length(something(findfirst(":rc", readname), 0:-1)) > 0
            @test align.value[1].strand == false
         else
            @test align.value[1].strand == true
         end

         best_ind = argmax(scores)

         ex_num = length(split(readname, '-', keepempty=false))
         @test length(align.value[best_ind].path) == ex_num

         count!( gquant, align.value[best_ind], 1.0 )

         curgene   = align.value[best_ind].path[1].gene
         firstnode = align.value[best_ind].path[1].node
         lastnode  = align.value[best_ind].path[end].node
         curgraph  = lib.graphs[ curgene ]

         cigar,positions = split(readname, '%', keepempty=false)[1:2]
         first,last   = split(positions, ',', keepempty=false) |> y->map(x->parse(Int,x), y)
         # Test SAM Format
         @test cigar == cigar_string( align.value[best_ind], curgraph, true, length(r.sequence) )[1]
         test_cigar,endpos = cigar_string( align.value[best_ind], curgraph, true, length(r.sequence) )
         #println(stderr, "cigar = $test_cigar")
         @test first == Int(curgraph.nodecoord[firstnode] + (align.value[best_ind].offset - curgraph.nodeoffset[firstnode]))
         #@test last  == Int(curgraph.nodecoord[lastnode] + (endpos - curgraph.nodeoffset[lastnode]))
         # test readlength = number of M and S entries in cigar
         # test SAM offset is correct for both '+' and '-' genes.
         # test that cigar reversal works for '-' strand genes.
      end

      @testset "SGAlignContainer Hashing" begin

         @test SGAlignNode(1,1,zero(SGAlignScore)) == SGAlignNode(1,1,one(SGAlignScore))
         @test hash( SGAlignNode(1,1,zero(SGAlignScore)) ) == hash( SGAlignNode(1,1,one(SGAlignScore)) )

         @test SGAlignNode[SGAlignNode(1,1,zero(SGAlignScore))] == SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))]
         @test hash( SGAlignNode[SGAlignNode(1,1,zero(SGAlignScore))] ) == hash( SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))] )

         @test SGAlignSingle(SGAlignNode[SGAlignNode(1,1,zero(SGAlignScore))]) == SGAlignSingle(SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))])
         @test hash( SGAlignSingle(SGAlignNode[SGAlignNode(1,1,zero(SGAlignScore))]) ) == hash( SGAlignSingle(SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))]) )

         a = map( x->SGAlignSingle(SGAlignNode[SGAlignNode(x,x,zero(SGAlignScore))]), collect(1:5))
         ab = map( x->SGAlignSingle(SGAlignNode[SGAlignNode(x,x,one(SGAlignScore))]), collect(1:5))
         b = map( x->SGAlignPaired(SGAlignNode[SGAlignNode(x,x,zero(SGAlignScore))], SGAlignNode[SGAlignNode(x,x,zero(SGAlignScore))]), collect(1:5))
         bb = map( x->SGAlignPaired(SGAlignNode[SGAlignNode(x,x,one(SGAlignScore))], SGAlignNode[SGAlignNode(x,x,one(SGAlignScore))]), collect(1:5))

         da = Dict{SGAlignSingle,Int}()
         db = Dict{SGAlignPaired,Int}()

         function settwice!( d, a, ab )
            for i in a
               d[i] = 1
            end
            for i in ab
               d[i] = 2
            end
            for k in keys(d)
               @test d[k] == 2
            end
         end
         settwice!( da, a, ab )
         settwice!( db, b, bb )
         @test length(keys(da)) == length(a)
         @test length(keys(db)) == length(b)

      end

      @testset "SGAlignPath Overlap" begin

         single     = SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))]
         single_two = SGAlignNode[SGAlignNode(1,2,one(SGAlignScore))]
         two_three  = SGAlignNode[SGAlignNode(1,2,one(SGAlignScore)), SGAlignNode(1,3,one(SGAlignScore))]
         larger     = SGAlignNode[SGAlignNode(1,1,one(SGAlignScore)), SGAlignNode(1,2,one(SGAlignScore)), SGAlignNode(1,3,one(SGAlignScore))]

         println(stderr, @which single in single)

         @test single in single
         @test !(single in single_two)
         @test single_two in two_three
         @test !(single in two_three)
         @test single in larger
         @test single_two in larger
         @test two_three in larger

      end

      @testset "Isoform Equivalence Classes" begin

         is = BitSet([1, 2, 5, 6])
         @test in( 1, 2, is ) == true
         @test in( 1, 3, is ) == false
         @test in( 2, 3, is ) == false
         @test in( 2, 4, is ) == false
         @test in( 2, 5, is ) == true
         @test in( 5, 6, is ) == true
         @test in( 1, 6, is ) == false
         for i in [1,2,5,6]
            @test (i in is) == true
         end
         path = SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000007, SGAlignScore(0x0a, 0x00, 0.0))]
         @test (path in is) == false
         path = SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000002, SGAlignScore(0x0a, 0x00, 0.0))]
         @test (path in is) == true
         path = SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000002, SGAlignScore(0x0a, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000003, SGAlignScore(0x0a, 0x00, 0.0))]
         @test (path in is) == false
         path = SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000002, SGAlignScore(0x0a, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000005, SGAlignScore(0x0a, 0x00, 0.0))]
         @test (path in is) == true

         @test length(gquant.tpm) == 7
         @test length(gquant.count) == 7
         @test length(gquant.length) == 7
         @test gquant.geneidx == [0, 2, 5, 6]
         build_equivalence_classes!( gquant, lib, assign_long=true )
         println(stderr, gquant)
         println(stderr, gquant.classes)

      end

      @testset "MultiMapping Equivalence Classes" begin

         multi  = MultiMapping{SGAlignSingle,DefaultCounter}()
         aligns = SGAlignment[SGAlignment(0x0000000e, 0x01, SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000007, SGAlignScore(0x0a, 0x00, 0.0))], false, true),
                              SGAlignment(0x00000001, 0x01, SGAlignNode[SGAlignNode(0x00000004, 0x00000001, SGAlignScore(0x02, 0x00, 0.0))], false, true)]

         push!( multi, aligns, 1.0, gquant, lib )
         println(stderr, multi.map)
         compats = collect(values(multi.map))
         @test compats[1].isdone == true
         @test sum(multi.iset) == 0

         multi  = MultiMapping{SGAlignSingle,DefaultCounter}()

         aligns = SGAlignment[SGAlignment(0x0000000e, 0x01, SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000002, SGAlignScore(0x0a, 0x00, 0.0))], false, true),
                              SGAlignment(0x00000001, 0x01, SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0))], false, true)]

         push!( multi, aligns, 1.0, gquant, lib )
         println(stderr, multi.map)
         compats = collect(values(multi.map))
         @test compats[1].isdone == false
         @test [3,4,5] == compats[1].class
         @test sum(multi.iset) == 0

         # Assignment
         compats[1].count = 10.0
         prev_node = get(gquant.quant[2].node[1])
         println(stderr, "prev_quant = $(gquant.quant[2])")
         interv = IntervalTrees.Interval{ExonInt}( 1, 2 )
         prev_edge = get(get( gquant.quant[2].edge, interv, IntervalValue(0,0,zero(DefaultCounter)) ).value)

         println(gquant.quant[2].edge)

         assign_ambig!( gquant, lib, multi )
         println(gquant.quant[2].edge)

         println(stderr, "cur_quant = $(gquant.quant[2])")
         @test get(gquant.quant[2].node[1]) == prev_node + 7.5
         @test get(get(gquant.quant[2].edge, interv, IntervalValue(0,0,default(DefaultCounter)) ).value) == prev_edge + 2.5
      end

      @testset "TPM" begin
         multi  = MultiMapping{SGAlignSingle,DefaultCounter}()
         calculate_tpm!( gquant, readlen=20 )
         set_gene_tpm!( gquant, lib )
         orig = deepcopy(gquant.tpm)
         iter = gene_em!( gquant, multi, sig=1, readlen=20, maxit=10000 )
         @test orig != gquant.tpm
         @test iter < 5000
         @test sum(gquant.tpm) == 1000000
      end

      @testset "Output Files" begin
         output_junctions( "test.jnc.gz", lib, gquant )
         output_stats( "test.map.gz", lib, gquant, param, "test_index.jls", 10, 10, 0, 20, "vTEST" )
         output_exons( "test_index.exons.gz", lib )
      end

      function parse_edge( str::S ) where S <: AbstractString
         s = split( str, ['-',':'] )
         (parse(Int, String(s[1])), parse(Int, String(s[2])), parse(Float64, String(s[3])))
      end

      @testset "Graph Enumeration" begin
         #= Real Test based on node 23 (here annotated as 1) in ENSMUSG00000038685
         inclusion edges:
         1-2, 2-3, 3-4, 4-5, 5-6, 6-7, 7-8, 8-9, 9-10
         exclusion edges:
         1-3, 2-4, 3-6, 3-10, 4-6, 8-10
         expected paths:
         =#
         expected_paths = [BitSet([1,2,3,4,5,6,7,8,9,10]) # (full)
         BitSet([1,3,4,5,6,7,8,9,10]) # (1-3)
         BitSet([1,2,4,5,6,7,8,9,10]) # (2-4)
         BitSet([1,2,3,6,7,8,9,10])   # (3-6)
         BitSet([1,3,6,7,8,9,10])     # (3-6)
         BitSet([1,2,3,10])           # (3-10)
         BitSet([1,3,10])             # (3-10)
         BitSet([1,2,3,4,6,7,8,9,10]) # (4-6)
         BitSet([1,3,4,6,7,8,9,10])   # (4-6)
         BitSet([1,2,4,6,7,8,9,10])   # (4-6)
         BitSet([1,2,3,4,5,6,7,8,10]) # (8-10)
         BitSet([1,3,4,5,6,7,8,10])   # (8-10)
         BitSet([1,2,4,5,6,7,8,10])   # (8-10)
         BitSet([1,2,3,6,7,8,10])     # (8-10)
         BitSet([1,3,6,7,8,10])       # (8-10)
         BitSet([1,2,3,4,6,7,8,10])   # (8-10)
         BitSet([1,3,4,6,7,8,10])     # (8-10)
         BitSet([1,2,4,6,7,8,10])]    # (8-10)
         edgestr = "1-2:616.0,1-3:4.0,2-3:569.0,2-4:20.0,3-4:629.0,3-6:1.0,3-10:3.0,4-5:1.0,4-6:664.0,5-6:5.0,6-7:789.0,7-8:790.0,8-9:2.0,8-10:606.0,9-10:15.0"
         edgespl  = split( edgestr, ',' )
         edges    = IntervalMap{Int,Float64}()
         psigraph = PsiGraph( Vector{Float64}(), Vector{Float64}(),
                              Vector{Float64}(), Vector{BitSet}(), 1, 10 )
         for s in split( edgestr, ',' )
            f,l,v = parse_edge( s )
            edges[(f,l)] = v
            push!( psigraph.psi, 0.0 )
            push!( psigraph.length, 1.0 )
            push!( psigraph.count, v )
            push!( psigraph.nodes, BitSet([f,l]) )
         end

         res = reduce_graph( psigraph )
         for p in expected_paths
            @test p in res.nodes
         end
         @test length(expected_paths) == length(res.nodes)
      end

      @testset "Event Building" begin
         out = IOBuffer(read=true,write=true)
         bs  = BufferedOutputStream(out)
         output_psi_header( bs )
         for g in 1:length(lib.graphs)
            name = lib.names[g]
            chr  = lib.info[g].name
            strand = lib.info[g].strand ? '+' : '-'
            _process_events( bs, lib.graphs[g], gquant.quant[g], (name,chr,strand), readlen=20, isnodeok=false )
         end
         flush(bs)
         seek(out,0)
         for l in eachline(out)
            spl = split( l, '\t' )
            @test length(spl) == 14
            println(stderr,l)
         end
      end
   end
end
