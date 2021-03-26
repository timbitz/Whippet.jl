# index.jl - tim.sterne.weiler@utoronto.ca, 1/28/16

abstract type SeqLibrary end

struct GraphLib <: SeqLibrary
   offset::Vector{CoordInt}
   names::Vector{GeneName}
   info::Vector{GeneInfo}
   lengths::Vector{Float64}
   graphs::Vector{SpliceGraph}
   edges::Edges
   index::FMIndex
   sorted::Bool
   kmer::Int64
end


function build_chrom_dict( ref::RefSet )
   ret = Dict{SeqName,Vector{GeneName}}() # refgenomeseq->geneid[]
   for g in keys(ref)
      chrom = ref[g].info.name
      if !haskey(ret, chrom)
         ret[chrom] = Vector{GeneName}()
      end
      push!(ret[chrom], g)
   end
   ret
end

function twobit_enc(seq)
   len = length(seq)
   ret = IntVector{2,UInt8}(len)
   for i in 1:len
      val = isambiguous(seq[i]) ? 0x00 : convert(UInt8, trailing_zeros(seq[i]))
      ret[i] = val
   end
   ret
end


function load_fasta( fhIter; verbose=false )
   seq = SGSequence(mutable=false)
   offset = Int[]
   names = SeqName[]
   for r in fhIter
      immutable!(r.seq)
      push!(names, r.name)
      push!(offset, length(seq)+1)
      println( stderr, r )
      @time seq *= r.seq
   end
   seq, offset, names
end

function single_genome_index!( fhIter; verbose=false )
   seq,offset,names = load_fasta( fhIter )
   @time fm = FMIndex(twobit_enc(seq), 4, r=4, program=:SuffixArrays, mmap=true)
   println( stderr, "Finished building Index..." )

   SeqLib(seq, offset, names, fm, false)
end


function trans_index!( fhIter, ref::RefSet; kmer=9 )
   seqdic  = build_chrom_dict( ref )
   xcript  = dna""
   xoffset = Vector{UInt64}()
   xgenes  = Vector{GeneName}()
   xinfo   = Vector{GeneInfo}()
   xlength = Vector{Float64}()
   xgraph  = Vector{SpliceGraph}()
   # set up splice graphs
   runoffset = 0
   num = 0
   for r in fhIter
      name = String(r.data[r.identifier])
      haskey(seqdic, name) || continue
      #BioSequences.immutable!(r.seq)
      sg = SGSequence( r.data[r.sequence] )

      println(stderr, "Building Splice Graphs for $( name ).." )
      @timer for g in seqdic[name]
         #println( stderr, "Building $g Splice Graph..." )
         curgraph = SpliceGraph( ref[g], sg, kmer )
         xcript  *= curgraph.seq
         push!(xgraph, curgraph)
         push!(xgenes, g)
         push!(xinfo, ref[g].info )
         push!(xlength, ref[g].length )
         push!(xoffset, runoffset)
         runoffset += length(curgraph.seq)
         num += 1
      end
   end
   if num > 1
      println( stderr, "Building full sg-index from $num genes..." )
      @time fm = FMIndex(twobit_enc(xcript), 4, r=1, program=:SuffixArrays, mmap=true)
   else
      error("ERROR: No genes in the GTF file matched chromosome names in the FASTA file!")
   end

   # clean up
   xcript = SGSequence()
   GC.gc()

   println( stderr, "Building edges.." )
   @time edges = build_edges( xgraph, kmer ) # TODO make variable kmer

   GraphLib( xoffset, xgenes, xinfo, xlength, xgraph, edges, fm, true, kmer )
end

function fasta_to_index( filename::String, ref::RefSet; kmer=9 )
   if isgzipped( filename )
      println(stderr, "Decompressing and Indexing $filename...")
      to_open = open( filename ) |> x->ZlibInflateInputStream(x, reset_on_end=true)
   else
      println(stderr, "Indexing $filename...")
      to_open = open( filename ) |> BufferedInputStream
   end
   # iterate through fasta entries
   index = @time trans_index!(FASTA.Reader( to_open ), ref, kmer=kmer)
   index
end
