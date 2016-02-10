using DataStructures
using IntervalTrees

#requires
include("types.jl")

#Base.convert{T<:AbstractString}(::Type{T}, seq::NucleotideSequence) = Base.convert(ASCIIString, [convert(Char, x) for x in seq])

typealias Refseqid    ASCIIString
typealias Txinfo      Tuple{Genename,Coordint,Coordint,Coordint} 
                       #   {Genename, TxStart, TxEnd, ExonCount}
typealias Geneinfo    Tuple{ASCIIString, Char}
                       #    Chrom/seqname, Strand '+'/'-'
typealias Coordtree IntervalTree{Coordint,Interval{Coordint}}

# Single Reftx entry
immutable Reftx
   info::Txinfo
   don::Coordtuple
   acc::Coordtuple
end

# Single Refgene entry
immutable Refgene
   info::Geneinfo
   #meanlen::Coordint
   don::Coordtuple # make Coordarray instead?
   acc::Coordtuple # TODO?
   txst::Coordtuple
   txen::Coordtuple
   exons::Coordtree
end

# Full Annotation Set
type Refset
   txset::Dict{Refseqid,Reftx}
   geneset::Dict{Genename,Refgene}
   genetotx::Dict{Genename,Vector{Refseqid}}
end


function Base.show(io::Base.IO, ref::Refset)
   Base.show(io, ref.txset)
   Base.show(io, ref.geneset)
   Base.show(io, ref.genetotx)
end


# this guy maps an array of strings or substrings to Int parsing and then into a tuple
# which is then sliced based on the l (left) and r (right) adjusters, c is the mathmatical adjuster to the data
parseSplice( subArray; l=0, r=0, c=0 ) = tuple(map( x->convert(Coordint,parse(Int,x)+c), subArray )...)[(1+l):(end-r)]

function unique_tuple{T}( tup1::Tuple{Vararg{T}}, tup2::Tuple{Vararg{T}})
   uniq = SortedSet{T,Base.Order.ForwardOrdering}()
   for i in tup1, j in tup2
     push!(uniq, i)
     push!(uniq, j)
   end
   tuple(uniq...)
end


# Load refflat file from filehandle
# Refflat format must be as expected from output of gtfToGenePred -ext
# Options:
#         if txbool=false, then return Refset with empty txset variable
function load_refflat( fh; txbool=false )
   
   # Temporary variables   
   gninfo = Dict{Genename,Geneinfo}()
   gndon = Dict{Genename,Coordtuple}()
   gnacc = Dict{Genename,Coordtuple}()
   gntxst = Dict{Genename,Coordtuple}()
   gntxen = Dict{Genename,Coordtuple}()
   gnlens = Dict{Genename,Coordtuple}()
   gnexons = Dict{Genename,Coordtree}()

   # Refset variables
   txset    = Dict{Refseqid,Reftx}()
   geneset  = Dict{Genename,Refgene}()
   genetotx = Dict{Genename,Vector{Refseqid}}()


   txnum = 75000
   sizehint!(txset, txnum)
   sizehint!(geneset, txnum >> 1)
   sizehint!(genetotx, txnum)

   tuppar( i; c=0 ) = tuple(convert(Coordint, parse(Int, i)+c))

   for l::ASCIIString in eachline(fh)
      (refid,chrom,strand, # NM_001177644,chr3,+
       txS,txE, # Int,Int
       cdS,cdE, # Int,Int
       exCnt, # Int,
       accCom,donCom, # 'Int,Int,Int','Int,Int,Int'
       _,gene) = split(l, '\t')

      exCnt = parse(UInt16, exCnt)

      if exCnt <= 2 continue end # no alternative splicing possible  TODO

      # get donor and acceptor splice sites, adjust for 0-based coords 
      don = split(donCom, '\,', keep=false) |> s->parseSplice(s, r=0)
      acc = split(accCom, '\,', keep=false) |> s->parseSplice(s, l=0, c=1)

      # Add original exons to interval tree-->
      for i in 1:length(don)
         if haskey(gnexons, gene)
            insval = Interval{Coordint}(acc[i],don[i])
            # make sure we are adding a unique value
            if !haskey(gnexons[gene], (acc[i],don[i]))
               push!(gnexons[gene], insval)
            end
         else
            insval = Interval{Coordint}(acc[i],don[i])
            gnexons[gene] = Coordtree()
            push!(gnexons[gene], insval)
         end
      end

      don = don[1:(end-1)] #ignore txStart and end
      acc = acc[2:end] # TODO make safe for single exon genes

      # set values
      txinfo = (gene,parse(Coordint, txS),
                     parse(Coordint, txE)-1,
                     exCnt)

      if txbool
         txset[refid] = Reftx( txinfo, don, acc )      
      end

      if haskey(genetotx, gene)
         chrom == gninfo[gene][1] || continue # can have only one chrom
         push!(genetotx[gene], refid)
         gndon[gene] = unique_tuple(gndon[gene], don)
         gnacc[gene] = unique_tuple(gnacc[gene], acc)
         gntxst[gene] = unique_tuple(gntxst[gene], tuppar(txS, c=1))
         gntxen[gene] = unique_tuple(gntxen[gene], tuppar(txE))
         #gnlens[gene] += genelen(don, acc)  # TODO finish
      else
         genetotx[gene] = Refseqid[refid]
         gndon[gene] = don
         gnacc[gene] = acc
         gninfo[gene] = (chrom,strand[1])
         gntxst[gene] = tuppar(txS, c=1)
         gntxen[gene] = tuppar(txE)
         #gnlens[gene] = genelen(don, acc)
      end
   end
   # now make Refset and add genes.
   for gene in keys(genetotx)
      geneset[gene] = Refgene( gninfo[gene], gndon[gene], gnacc[gene],
                               gntxst[gene], gntxen[gene], gnexons[gene] )
   end

   if !txbool
      genetotx = Dict{Genename,Vector{Refseqid}}()
      gc()
   end

   return Refset( txset, geneset, genetotx )
end


function main()
   
   println(STDERR, "Loading Refflat file...")
   fh = open("$(pwd())/../genome/genes.flat", "r")
   @time ref = load_refflat(fh)
   close(fh)

   println(STDERR, "Saving gene annotations...")
   open("$(pwd())/../index/reflat.jls", "w+") do fh
      @time serialize(fh, ref)
   end 
   println(STDERR, "Loading annotation index...")
   @time ref = open(deserialize, "$(pwd())/../index/reflat.jls")
   #=
   Loading Refflat file...
    45.661196 seconds (149.57 M allocations: 3.433 GB, 4.67% gc time)
   Saving gene annotations...
     2.310959 seconds (3.66 M allocations: 126.598 MB, 2.11% gc time)
   Loading annotation index...
     1.510497 seconds (6.74 M allocations: 175.154 MB, 18.79% gc time)
   =#
   gc()
   return ref  

   println(STDERR, "Loading genome index...")
   @time genome = open(deserialize, "$(pwd())/index/genome.jls")
   @time junclib_d,junclib_a = donor_junc_table( genome, ref, 9 )
   open("$(pwd())/index/donor.jls", "w+") do fh
      @time serialize(fh, junclib_d)
      @time serialize(fh, junclib_a)
   end
   gc()
   sleep(10)
end

refflat = main()
