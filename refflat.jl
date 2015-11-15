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
typealias Geneinfo    Tuple{Genename,ASCIIString,Char,Coordint,Coordint,Coordint} 
                        # {Genename, Chrom, Strand, TxStart, TxEnd, ExonCount}
typealias Splicetuple Tuple{Vararg{Coordint}}

immutable Refflat
   txinfo::Dict{Refseqid,Geneinfo}
   txdon::Dict{Refseqid,Splicetuple}
   txacc::Dict{Refseqid,Splicetuple}
   idtogene::Dict{Genename,Vector{Refseqid}}
   gndon::Dict{Genename,Splicetuple}
   gnacc::Dict{Genename,Splicetuple}
end

function Base.show(io::IO, ref::Refflat)
   Base.show(io, ref.txinfo)
   Base.show(io, ref.txdon)
   Base.show(io, ref.txacc)
   Base.show(io, ref.idtogene)
end

parseSplice( subArray ) = tuple(map( x->parse(Coordint,x), subArray )...)

function unique_tuple{T}( tup1::Tuple{Vararg{T}}, tup2::Tuple{Vararg{T}})
   uniq = SortedSet{T,Base.Order.ForwardOrdering}()
   for i in tup1, j in tup2
     push!(uniq, i)
     push!(uniq, j)
   end
   tuple(uniq...)
end

function load_refflat( fh )
   txinfo = Dict{Refseqid,Geneinfo}()
   txdon = Dict{Refseqid,Splicetuple}()
   txacc = Dict{Refseqid,Splicetuple}()
   genetotx = Dict{Genename,Vector{Refseqid}}()

   gndon = Dict{Genename,Splicetuple}()
   gnacc = Dict{Genename,Splicetuple}()

   txnum = 75000
   sizehint!(txinfo, txnum)
   sizehint!(txdon, txnum)
   sizehint!(txacc, txnum)
   sizehint!(genetotx, txnum)

   for l::ASCIIString in eachline(fh)
      (refid,chrom,strand, # NM_001177644,chr3,+
       txS,txE, # Int,Int
       cdS,cdE, # Int,Int
       exCnt, # Int,
       accCom,donCom, # 'Int,Int,Int','Int,Int,Int'
       _,gene) = split(l, '\t')
      don = split(donCom, '\,', keep=false) |> parseSplice
      acc = split(accCom, '\,', keep=false) |> parseSplice

      # set values
      txinfo[refid] = (gene,chrom,strand[1],
                                  parse(Coordint, txS),
                                  parse(Coordint, txE),
                                  parse(Coordint,exCnt))
      txdon[refid] = don
      txacc[refid] = acc
      if haskey(genetotx, gene)
         push!(genetotx[gene], refid)
         gndon[gene] = unique_tuple(gndon[gene], don)
         gnacc[gene] = unique_tuple(gnacc[gene], acc)
      else
         genetotx[gene] = Refseqid[refid]
         gndon[gene] = don
         gnacc[gene] = acc
      end
   end
   return Refflat( txinfo, txdon, txacc, genetotx, gndon, gnacc )
end

immutable Acceptor{K}
#   kmer::DNAKmer{K}
   microexon::UInt8
   gene::Genename
end

function doner_junc_table( genome, ref, k::Integer )
    doners = Vector{Vector{Acceptor{k}}}(4^k)
    for gene in keys(ref.genetotx)
       for tx in ref.genetotx[gene]
                    
       end
    end   
end

function main()
   println(STDERR, "Loading Refflat file...")
   fh = open("$(pwd())/genome/genes.flat", "r")
   @time ref = load_refflat(fh)
   close(fh)

   println(STDERR, "Saving gene annotations...")
   open("$(pwd())/index/reflat.jls", "w+") do fh
      @time serialize(fh, ref)
   end

   println(STDERR, "Loading annotation index...")
   @time trans = open(deserialize, "$(pwd())/index/reflat.jls")
   println(STDERR, "Loading genome index...")
   @time genome = open(deserialize, "$(pwd())/index/genome.jls")
   sleep(10)
end

main()
