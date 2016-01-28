using GZip
#using Bio.Seq
using FMIndexes
using DataStructures

#Base.convert{T<:AbstractString}(::Type{T}, seq::NucleotideSequence) = Base.convert(ASCIIString, [convert(Char, x) for x in seq])

const Mb = 1_000_000
const GENOMESIZE = 3235Mb

# ALIAS NEW TYPES FOR INCREASED CODE READABILITY
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
typealias Coordarray Vector{Coordint}
typealias Microsize  UInt8
typealias Microtuple Tuple{Vararg{Microsize}}

# Single Reftx entry
immutable Reftx
   info::Txinfo
   don::Coordtuple
   acc::Coordtuple
end

# Single Refgene entry
immutable Refgene
   info::Geneinfo
   meanlen::Coordint
   don::Coordtuple # make Coordarray instead?
   acc::Coordtuple # TODO?
   txst::Coordtuple
   txen::Coordtuple
   dondic::Dict{Coordint,Coordtuple}
end

# Full Annotation Set
type Refset
   txset::Dict{Refseqid,Reftx}
   geneset::Dict{Genename,Refgene}
   genetotx::Dict{Genename,Vector{Refseqid}}
end

"""
function Base.show(io::IO, ref::Refflat)
   Base.show(io, ref.txinfo)
   Base.show(io, ref.txdon)
   Base.show(io, ref.txacc)
   Base.show(io, ref.genetotx)
end
"""

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

# Load refflat file from filehandle
# Refflat format must be as expected from output of gtfToGenePred -ext
function load_refflat( fh )
   
   gninfo = Dict{Genename,Geneinfo}()
   gndon = Dict{Genename,Coordtuple}()
   gnacc = Dict{Genename,Coordtuple}()
   gntxst = Dict{Genename,Coordtuple}()
   gntxen = Dict{Genename,Coordtuple}()
#   gnmic = Dict{Genename,Microtuple}()
   gnlens = Dict{Genename,Coordtuple}()

   # Refset variables
   txset    = Dict{Refseqid,Reftx}()
   geneset  = Dict{Genename,Refgene}()
   genetotx = Dict{Genename,Vector{Refseqid}}()


   txnum = 75000
   sizehint!(txset, txnum)
   sizehint!(geneset, txnum >> 1)
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

      # get donor and acceptor splice sites, ignore txStart/End (l/r), 
      #                                      and adjust for 0-based coords
      don = split(donCom, '\,', keep=false) |> s->parseSplice(s, r=1, c=-1)
      acc = split(accCom, '\,', keep=false) |> s->parseSplice(s, l=1)

      # set values
      txinfo = (gene,parse(Coordint, txS),
                     parse(Coordint, txE),
                     exCnt)

      txset[refid] = Reftx( txinfo, don, acc )      

      if haskey(genetotx, gene)
         push!(genetotx[gene], refid)
         gndon[gene] = unique_tuple(gndon[gene], don)
         gnacc[gene] = unique_tuple(gnacc[gene], acc)
         gntxst[gene] = unique_tuple(gntxst[gene], tuppar(txS))
         gntxen[gene] = unique_tuple(gntxen[gene], tuppar(txE))
         gnlens[gene] += genelen(don, acc)  #FINISH
      else
         genetotx[gene] = Refseqid[refid]
         gndon[gene] = don
         gnacc[gene] = acc
         gninfo[gene] = (chrom,strand[1])
         gntxst[gene] = tuppar(txS)
         gntxen[gene] = tuppar(txE)
         gnlens[gene] = genelen(don, acc)
      end
   end
   # now make Refset and add genes.
   for gene in keys(genetotx)
      dondic = Dict{Coordint,Coordtuple}()
      # make dondic
      geneset[gene] = Refgene( gninfo[gene], gndon[gene], gnacc[gene],
                               gntxst[gene], gntxen[gene], dondic )
   end

   return Refset( txset, geneset, genetotx )
end


#### Splice k-mer graphs.
immutable Acceptor{K}
   kmer::DNAKmer{K}
   microexon::UInt8
   gene::Genename
end

function setacceptor!{K}( donorvec::Vector, left::DNAKmer{K}, right::DNAKmer{K}, gene::Genename )
   lind = Int(UInt64(left)) + 1
   if !isdefined( donorvec, lind )
      donorvec[lind] = Vector{Acceptor{K}}()
   end
   push!(donorvec[lind], Acceptor{K}(right, UInt8(0), gene) )
end

function push_or_set!{T}( array::Vector, index::DNAKmer, value::T)
   ind = Int(UInt64(index)) + 1
   if !isdefined( array, ind )
      array[ind] = Set{T}()
   end
   push!( array[ind], value )
end

function Base.intersect{T}( tup1::Tuple{Vararg{T}}, tup2::Tuple{Vararg{T}} )
   res = Vector{T}()
   for a in tup1, b in tup2
      if a == b
         push!(res, a)
      end
   end
   res
end

# Intersect two sorted arrays.
function intersect_sorted{T}( arrA::Vector{T}, arrB::Vector{T} )
   res = Vector{T}()
   i,j = 1,1
   while i <= length(arrA) && j <= length(arrB)
      if arrA[i] > arrB[j]
         j += 1
      elseif arrA[i] < arrB[j]
         i += 1
      else
         push!( res, arrA[i] )
         i += 1
         j += 1
      end
   end
   res
end 

# Binary search.
function search_sorted{T}( arr::Vector{T}, elem::T, low=1, high=length(arr)+1 )
   low == high && return(-1)
   mid = ((high - low) >> 1) + low
   arr[mid] == elem && return(mid)
   if arr[mid] > elem
      ret = search_sorted(arr, elem, low, mid)
   else
      ret = search_sorted(arr, elem, mid+1, high)
   end
   ret
end

# build k-mer arrays for splice junctions
function donor_junc_table( genome, ref, k::Integer )
   donors = Vector{Vector{Genename}}(4^k)
   accept = Vector{Vector{Genename}}(4^k)

   for gene in keys(ref.genetotx)
      chrom,strand = ref.gninfo[gene]
      offset = getoffset( genome, chrom )
      (offset == -1) && continue # chrom not present

      function build_kmer_gene( coords::Tuple, kmers::Vector; ladj=0, radj=0 )
         for c in coords
            coor = c+offset
            try
               curkmer = DNAKmer{k}(genome.seq[(coor+ladj):(coor+radj)])
               ind = Int(UInt64(curkmer))+1
               if !isdefined(kmers, ind)
                  kmers[ind] = Vector{typeof(gene)}()
               end
               push!(kmers[ind], gene)
               kmers[ind] = sort(kmers[ind]) |> unique
            catch e
               #println(e) # these should just be N containing
               continue
            end
         end
      end
      build_kmer_gene( ref.gndon[gene], donors, ladj=(1-k) )
      build_kmer_gene( ref.gnacc[gene], accept, radj=(k-1) )

   end
   donors,accept
end
      #=
      for d in ref.gndon[gene]#, a in ref.gnacc[gene]
         don = d+offset#,a+offset
         try
            left = DNAKmer{k}(genome.seq[(don-k+1):don])
            #right = DNAKmer{k}(genome.seq[acc:(acc+k-1)])
            #push_or_set!( donors, left, gene )
            lind = Int(UInt64(left))+1
            donors[lind] = isdefined(donors, lind) ? unique_tuple(donors[lind], tuple(gene)) : tuple(gene)
         catch e
            println(e)
            continue
         end
         #println("$gene $chrom $strand $d $k $curk")
         #print(">$chrom-$don-$acc\n$(convert(AbstractString,left)*convert(AbstractString,right))\n")
      end
      =#

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

   include("index.jl")
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
