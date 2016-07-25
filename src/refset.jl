
# Single RefTx entry
immutable RefTx
   info::TxInfo
   don::CoordTuple
   acc::CoordTuple
   length::Float64
end

# Single RefGene entry
immutable RefGene
   info::GeneInfo
   don::CoordTuple 
   acc::CoordTuple 
   txst::CoordTuple
   txen::CoordTuple
   orfst::CoordTuple
   orfen::CoordTuple
   exons::CoordTree
   length::Float64
end

# Full Annotation Set
type RefSet
   txset::Dict{RefSeqId,RefTx}
   geneset::Dict{GeneName,RefGene}
   genetotx::Dict{GeneName,Vector{RefSeqId}}
end


function Base.show(io::Base.IO, ref::RefSet)
   Base.show(io, ref.txset)
   Base.show(io, ref.geneset)
   Base.show(io, ref.genetotx)
end


# this guy maps an array of strings or substrings to Int parsing and then into a tuple
# which is then sliced based on the l (left) and r (right) adjusters, c is the mathmatical adjuster to the data
parse_splice( subarray; l=0, r=0, c=0 ) = tuple(map( x->convert(CoordInt,parse(Int,x)+c), subarray )...)[(1+l):(end-r)]
parse_coordint( i ) = convert(CoordInt, parse(Int, i))
parse_coordtup( i; c=0 ) = tuple(parse_coordint(i)+c)

unique_tuple( tup1::Tuple{}, tup2::Tuple{} ) = ()
unique_tuple{T}( tup1::Tuple{Vararg{T}}, tup2::Tuple{} ) = unique_tuple(tup1, tup1)
unique_tuple{T}( tup1::Tuple{}, tup2::Tuple{Vararg{T}} ) = unique_tuple(tup2, tup2)

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
#         if txbool=false, then return RefSet with empty txset variable
function load_refflat( fh; txbool=true )
   
   # Temporary variables   
   gninfo   = Dict{GeneName,GeneInfo}()
   gndon    = Dict{GeneName,CoordTuple}()
   gnacc    = Dict{GeneName,CoordTuple}()
   gntxst   = Dict{GeneName,CoordTuple}()
   gntxen   = Dict{GeneName,CoordTuple}()
   gnorfst  = Dict{GeneName,CoordTuple}()
   gnorfen  = Dict{GeneName,CoordTuple}()
   gnlens   = Dict{GeneName,CoordTuple}()
   gnexons  = Dict{GeneName,CoordTree}()
   gnlen    = Dict{GeneName,Float64}()
   gncnt    = Dict{GeneName,Int}()

   # RefSet variables
   txset    = Dict{RefSeqId,RefTx}()
   geneset  = Dict{GeneName,RefGene}()
   genetotx = Dict{GeneName,Vector{RefSeqId}}()


   txnum = 75000
   sizehint!(txset, txnum)
   sizehint!(geneset, txnum >> 1)
   sizehint!(genetotx, txnum)


   for l in eachline(fh)

      (l[1] == "#") && continue # ignore comment lines
      
      (refid,chrom,strand, # NM_001177644,chr3,+
       txS,txE, # Int,Int
       cdS,cdE, # Int,Int
       exCnt, # Int,
       accCom,donCom, # 'Int,Int,Int','Int,Int,Int'
       _,gene) = split(l, '\t')

      exCnt = parse(UInt16, exCnt)

      txlen = 0

      # get donor and acceptor splice sites, adjust for 0-based coords 
      don = split(donCom, '\,', keep=false) |> s->parse_splice(s, r=0)
      acc = split(accCom, '\,', keep=false) |> s->parse_splice(s, l=0, c=1)

      # Add original exons to interval tree-->
      for i in 1:length(don)
         insval = Interval{CoordInt}(acc[i],don[i])
         txlen += don[i] - acc[i] + 1
         if haskey(gnexons, gene)
            # make sure we are adding a unique value
            if !haskey(gnexons[gene], (acc[i],don[i]))
               push!(gnexons[gene], insval)
            end
         else
            gnexons[gene] = CoordTree()
            push!(gnexons[gene], insval)
         end
      end

      don = don[1:(end-1)] #ignore txStart and end
      acc = acc[2:end] 

      # set values
      txinfo = TxInfo(gene,parse(CoordInt, txS),
                     parse(CoordInt, txE)-1,
                     exCnt)

      if txbool
         txset[refid] = RefTx( txinfo, don, acc, txlen )      
      end

      if haskey(genetotx, gene)
         chrom == gninfo[gene].name || continue # can have only one chrom
         push!(genetotx[gene], refid)
         gndon[gene]  = unique_tuple(gndon[gene], don)
         gnacc[gene]  = unique_tuple(gnacc[gene], acc)
         gntxst[gene] = unique_tuple(gntxst[gene], parse_coordtup(txS, c=1))
         gntxen[gene] = unique_tuple(gntxen[gene], parse_coordtup(txE))
         gnlen[gene] += txlen
         gncnt[gene] += 1
      else
         genetotx[gene] = RefSeqId[refid]
         gndon[gene]  = don
         gnacc[gene]  = acc
         gninfo[gene] = GeneInfo(chrom,strand[1])
         gntxst[gene] = parse_coordtup(txS, c=1)
         gntxen[gene] = parse_coordtup(txE)
         gnlen[gene]  = txlen
         gncnt[gene]  = 1
      end
   end
   # now make RefSet and add genes.
   for gene in keys(genetotx)
      geneset[gene] = RefGene( gninfo[gene],  
                               gndon[gene],   gnacc[gene],
                               gntxst[gene],  gntxen[gene], 
                               gnexons[gene], 
                               gnlen[gene] /  gncnt[gene] )
   end

   if !txbool
      genetotx = Dict{GeneName,Vector{RefSeqId}}()
      gc()
   end

   return RefSet( txset, geneset, genetotx )
end

function fetch_meta( var::String, meta; off=1 )
   val = "_" # empty
   for i in off:length(meta)
      if meta[i] == var && i <= length(meta)
         return string(meta[i+1]),i+1
      end
   end
   val,0
end

function load_gtf( fh; txbool=true )

   # Temporary variables   
   gninfo   = Dict{GeneName,GeneInfo}()
   gndon    = Dict{GeneName,CoordTuple}()
   gnacc    = Dict{GeneName,CoordTuple}()
   gntxst   = Dict{GeneName,CoordTuple}()
   gntxen   = Dict{GeneName,CoordTuple}()
   gnorfst  = Dict{GeneName,CoordTuple}()
   gnorfen  = Dict{GeneName,CoordTuple}()
   gnlens   = Dict{GeneName,CoordTuple}()
   gnexons  = Dict{GeneName,CoordTree}()
   gnlen    = Dict{GeneName,Float64}()
   gncnt    = Dict{GeneName,Int}()

   # RefSet variables
   txset    = Dict{RefSeqId,RefTx}()
   geneset  = Dict{GeneName,RefGene}()
   genetotx = Dict{GeneName,Vector{RefSeqId}}()


   txnum = 75000
   sizehint!(txset, txnum)
   sizehint!(geneset, txnum >> 1)
   sizehint!(genetotx, txnum)

   curtran = ""
   curgene = ""

   trandon = Vector{CoordInt}()
   tranacc = Vector{CoordInt}()
   txlen   = 0

   for l in eachline(fh)

      (l[1] == "#") && continue # ignore comment lines

      (chrom, src, entrytype,
       st, en, _, 
       strand, _,
       meta) = split(l, '\t')

      (entrytype != "exon") || continue

      metaspl = split(meta, [' ',';','"'], keep=false)

      geneid,pos  = fetch_meta( "gene_id", metaspl )
      tranid,pos  = fetch_meta( "transcript_id", metaspl )
      genesym     = fetch_meta( "gene_name", metaspl )

      if tranid != curtran
         # clean up tran
         txinfo = TxInfo(curgene, minimum(tranacc),
                                  minimum(trandon),
                                  length(tranacc))

         don = tuple(trandon[1:(end-1)]...)
         acc = tuple(tranacc[2:end]...)

         if txbool
            txset[refid] = RefTx( txinfo, don, acc, txlen )
         end

         curtran = tranid

         if geneid != curgene
         # clean up gene

            curgene = geneid
         end
      end
 
      push!(tranacc, parse_coordint(st)) 
      push!(trandon, parse_coordint(en))
      txlen += trandon[end] - tranacc[end]

   end
end

