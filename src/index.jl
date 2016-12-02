# index.jl - tim.sterne.weiler@utoronto.ca, 1/28/16

# requires

abstract SeqLibrary

immutable GraphLib <: SeqLibrary
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


function build_offset_dict{I <: Integer, 
                           S <: AbstractString}( offset::Vector{I}, names::Vector{S} )
   ret = Dict{S,I}()
   for i in 1:length(names)
      cname,coffset = names[i],offset[i]
      ret[cname] = coffset
   end
   ret
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

# encode 3-bit sequence with L,R,S,N,A,T,G,C
function threebit_enc(seq)
   len = length(seq)
   ret = IntVector{3,UInt8}(len)
   for i in 1:len
      ret[i] = convert(UInt8, seq[i])
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
      println( STDERR, r )
      @time seq *= r.seq
   end
   seq, offset, names
end

function single_genome_index!( fhIter; verbose=false )
   seq,offset,names = load_fasta( fhIter )
   @time fm = FMIndex(twobit_enc(seq), 4, r=4, program=:SuffixArrays, mmap=true)
   println( STDERR, "Finished building Index..." )

   SeqLib(seq, offset, names, fm, false)
end


function trans_index!( fhIter, ref::RefSet; kmer=9 )
   seqdic  = build_chrom_dict( ref )
   xcript  = sg""
   xoffset = Vector{UInt64}()
   xgenes  = Vector{GeneName}()
   xinfo   = Vector{GeneInfo}()
   xlength = Vector{Float64}()
   xgraph  = Vector{SpliceGraph}()
   # set up splice graphs
   runoffset = 0
   for r in fhIter
      haskey(seqdic, r.name) || continue
      #Bio.Seq.immutable!(r.seq)
      sg = SGSequence( r.seq )

      r.seq = dna""
      gc() # free

      println(STDERR, "Building Splice Graphs for $( r.name ).." )
      for g in seqdic[r.name]
         #println( STDERR, "Building $g Splice Graph..." )
         curgraph = SpliceGraph( ref[g], sg, kmer )
         xcript  *= curgraph.seq
         push!(xgraph, curgraph)
         push!(xgenes, g)
         push!(xinfo, ref[g].info )
         push!(xlength, ref[g].length )
         push!(xoffset, runoffset)
         runoffset += length(curgraph.seq) 
      end
   end
   println( STDERR, "Building full sg-index.." )
   @time fm = FMIndex(threebit_enc(xcript), 8, r=1, program=:SuffixArrays, mmap=true) 

   # clean up
   xcript = sg""
   gc()

   println( STDERR, "Building edges.." ) 
   @time edges = build_edges( xgraph, kmer ) # TODO make variable kmer

   GraphLib( xoffset, xgenes, xinfo, xlength, xgraph, edges, fm, true, kmer )
end

function fasta_to_index( filename::String, ref::RefSet; kmer=9 )
   if isgzipped( filename )
      println(STDERR, "Decompressing and Indexing $filename...")
      to_open = open( filename ) |> x->ZlibInflateInputStream(x, reset_on_end=true)
   else
      println(STDERR, "Indexing $filename...")
      to_open = open( filename ) |> BufferedInputStream
   end
   # iterate through fasta entries
   index = @time trans_index!(FASTAReader{Bio.Seq.ReferenceSequence}( to_open, nothing ), ref, kmer=kmer)
   index
end

