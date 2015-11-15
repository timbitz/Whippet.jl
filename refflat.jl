using GZip
using Bio.Seq
using FMIndexes
using DataStructures

include("index.jl")

const Mb = 1000000
const GENOMESIZE = 3235Mb

if GENOMESIZE < typemax(UInt32) 
   typealias Coordint UInt32 
else 
   typealias Coordint UInt64
end
typealias Genename    ASCIIString
typealias Refseqid    ASCIIString
typealias Txinfo      Tuple{Genename,Coordint,Coordint,Coordint} 
                       #   {Genename, TxStart, TxEnd, ExonCount}
typealias Geneinfo    Tuple{ASCIIString, Char}
                       #    Chrom/seqname, Strand '+'/'-'
typealias Coordtuple Tuple{Vararg{Coordint}}

immutable Refflat
   txinfo::Dict{Refseqid,Txinfo}
   txdon::Dict{Refseqid,Coordtuple}
   txacc::Dict{Refseqid,Coordtuple} 
   genetotx::Dict{Genename,Vector{Refseqid}} # Array of Tx Refseqids
   gninfo::Dict{Genename,Geneinfo} # Seqname, Chrom
   gndon::Dict{Genename,Coordtuple} # Sorted tuple of exon ends
   gnacc::Dict{Genename,Coordtuple} # Sorted tuple of exon starts
   gntxst::Dict{Genename,Coordtuple} # Sorted tuple of tx starts
   gntxen::Dict{Genename,Coordtuple} # Sorted tuple of tx ends
end

function Base.show(io::IO, ref::Refflat)
   Base.show(io, ref.txinfo)
   Base.show(io, ref.txdon)
   Base.show(io, ref.txacc)
   Base.show(io, ref.idtogene)
end

parseSplice( subArray; l=0, r=0 ) = tuple(map( x->parse(Coordint,x), subArray )...)[(1+l):(end-r)]

function unique_tuple{T}( tup1::Tuple{Vararg{T}}, tup2::Tuple{Vararg{T}})
   uniq = SortedSet{T,Base.Order.ForwardOrdering}()
   for i in tup1, j in tup2
     push!(uniq, i)
     push!(uniq, j)
   end
   tuple(uniq...)
end

function find_microexons{T}( don::Tuple{Vararg{T}}, acc::Tuple{Vararg{T}})
   mics = Vector{Tuple{T,T,T}}()
   for i in 1:(length(acc)-1)
      exonlen = don[i+1] - acc[i]
      if exonlen < 16   #TODO make sure 1-based and 0-based coordinate math works!
         push!(mics, (acc[i], don[i+1], convert(T, i)))
      end
   end
   mics
end

function load_refflat( fh )
   txinfo = Dict{Refseqid,Txinfo}()
   txdon = Dict{Refseqid,Coordtuple}()
   txacc = Dict{Refseqid,Coordtuple}()
   genetotx = Dict{Genename,Vector{Refseqid}}()
   gninfo = Dict{Genename,Geneinfo}()
   gndon = Dict{Genename,Coordtuple}()
   gnacc = Dict{Genename,Coordtuple}()
   gntxst = Dict{Genename,Coordtuple}()
   gntxen = Dict{Genename,Coordtuple}()

   txnum = 75000
   sizehint!(txinfo, txnum)
   sizehint!(txdon, txnum)
   sizehint!(txacc, txnum)
   sizehint!(genetotx, txnum)

   tuppar( i ) = tuple(parse(Coordint, i))

   for l::ASCIIString in eachline(fh)
      (refid,chrom,strand, # NM_001177644,chr3,+
       txS,txE, # Int,Int
       cdS,cdE, # Int,Int
       exCnt, # Int,
       accCom,donCom, # 'Int,Int,Int','Int,Int,Int'
       _,gene) = split(l, '\t')

      exCnt = parse(UInt16, exCnt)

      if exCnt <= 2 continue end # no alternative splicing possible

      # get donor and acceptor splice sites, ignore txStart/End
      don = split(donCom, '\,', keep=false) |> s->parseSplice(s, r=1)
      acc = split(accCom, '\,', keep=false) |> s->parseSplice(s, l=1)

      # set values
      txinfo[refid] = (gene,parse(Coordint, txS),
                            parse(Coordint, txE),
                            exCnt)
      txdon[refid] = don
      txacc[refid] = acc

      if haskey(genetotx, gene)
         push!(genetotx[gene], refid)
         gndon[gene] = unique_tuple(gndon[gene], don)
         gnacc[gene] = unique_tuple(gnacc[gene], acc)
         gntxst[gene] = unique_tuple(gntxst[gene], tuppar(txS))
         gntxen[gene] = unique_tuple(gntxen[gene], tuppar(txE))
      else
         genetotx[gene] = Refseqid[refid]
         gndon[gene] = don
         gnacc[gene] = acc
         gninfo[gene] = (chrom,strand[1])
         gntxst[gene] = tuppar(txS)
         gntxen[gene] = tuppar(txE)
      end
   end
   return Refflat( txinfo, txdon, txacc, genetotx, gninfo, gndon, gnacc, gntxst, gntxen )
end


#### Splice k-mer graphs.
immutable Acceptor{K}
   kmer::DNAKmer{K}
   microexon::UInt8
   gene::Genename
end

function donor_junc_table( genome, ref, k::Integer )
   doners = Vector{Vector{Acceptor{k}}}(4^k)
   for gene in keys(ref.genetotx)
      chrom,strand = ref.gninfo[gene]
      offset = getoffset( genome, chrom )
      (offset == -1) && continue # chrom not present
      for d in ref.gndon[gene]
         coor = d+offset
         curk = DNAKmer{20}(genome.seq[(coor-20+1):coor])
         println("$gene $chrom $strand $d $k $curk")
      end
      break
   end   
end

function main()
"""   println(STDERR, "Loading Refflat file...")
   fh = open("$(pwd())/genome/genes.flat", "r")
   @time ref = load_refflat(fh)
   close(fh)

   println(STDERR, "Saving gene annotations...")
   open("$(pwd())/index/reflat.jls", "w+") do fh
      @time serialize(fh, ref)
   end
"""
   println(STDERR, "Loading annotation index...")
   @time ref = open(deserialize, "$(pwd())/index/reflat.jls")
   println(STDERR, "Loading genome index...")
   @time genome = open(deserialize, "$(pwd())/index/genome.jls")
   sleep(10)
   donor_junc_table( genome, ref, 8 )
end

main()
