
# Single RefTx entry
struct RefTx
   info::TxInfo
   don::CoordTuple
   acc::CoordTuple
   length::Float64
end

# Single RefGene entry
struct RefGene
   info::GeneInfo
   don::CoordTuple 
   acc::CoordTuple 
   txst::CoordTuple
   txen::CoordTuple
   exons::CoordTree
   length::Float64
   reftx::Vector{RefTx}
end

const RefSet = Dict{GeneName,RefGene}

# this guy maps an array of strings or substrings to Int parsing and then into a tuple
# which is then sliced based on the l (left) and r (right) adjusters, c is the mathmatical adjuster to the data
parse_splice( subarray; l=0, r=0, c=0 ) = tuple(map( x->convert(CoordInt,parse(Int,x)+c), subarray )...)[(1+l):(end-r)]
parse_coordint( i ) = convert(CoordInt, parse(Int, i))
parse_coordtup( i; c=0 ) = tuple(parse_coordint(i)+CoordInt(c))

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
   gnlens   = Dict{GeneName,CoordTuple}()
   gnexons  = Dict{GeneName,CoordTree}()
   gnreftx  = Dict{GeneName,Vector{RefTx}}()
   gnlen    = Dict{GeneName,Float64}()
   gncnt    = Dict{GeneName,Int}()

   # RefSet variables
   geneset  = Dict{GeneName,RefGene}()


   txnum = 75000
   sizehint!(geneset, txnum >> 1)

   for l in eachline(fh)

      (l[1] == '#') && continue # ignore comment lines
      
      (refid,chrom,strand, # NM_001177644,chr3,+
       txS,txE, # Int,Int
       cdS,cdE, # Int,Int
       exCnt, # Int,
       accCom,donCom, # 'Int,Int,Int','Int,Int,Int'
       _,gene) = split(chomp(l), '\t')

      exCnt = parse(UInt16, exCnt)

      txlen = 0

      # get donor and acceptor splice sites, adjust for 0-based coords 
      trandon = split(donCom, '\,', keep=false) |> s->parse_splice(s, r=0)
      tranacc = split(accCom, '\,', keep=false) |> s->parse_splice(s, l=0, c=1)

      # Add original exons to interval tree-->
      for i in 1:length(trandon)
         insval = Interval{CoordInt}(tranacc[i],trandon[i])
         txlen += trandon[i] - tranacc[i] + 1
         if haskey(gnexons, gene)
            # make sure we are adding a unique value
            if !haskey(gnexons[gene], (tranacc[i],trandon[i]))
               push!(gnexons[gene], insval)
            end
         else
            gnexons[gene] = CoordTree()
            push!(gnexons[gene], insval)
         end
      end

      don = trandon[1:(end-1)] #ignore txStart and end
      acc = tranacc[2:end] 

      # set values
      txinfo = TxInfo(refid,parse(CoordInt, txS)+1,
                            parse(CoordInt, txE),
                            exCnt)

      if haskey(gninfo, gene)
         chrom == gninfo[gene].name || continue # can have only one chrom
         gndon[gene]  = unique_tuple(gndon[gene], don)
         gnacc[gene]  = unique_tuple(gnacc[gene], acc)
         gntxst[gene] = unique_tuple(gntxst[gene], parse_coordtup(txS, c=1))
         gntxen[gene] = unique_tuple(gntxen[gene], parse_coordtup(txE))
         gnlen[gene] += txlen
         gncnt[gene] += 1
         push!( gnreftx[gene], RefTx( txinfo, trandon, tranacc, txlen ) )
      else
         gndon[gene]    = don
         gnacc[gene]    = acc
         gninfo[gene]   = GeneInfo(gene, chrom, strand[1])
         gntxst[gene]   = parse_coordtup(txS, c=1)
         gntxen[gene]   = parse_coordtup(txE)
         gnlen[gene]    = txlen
         gncnt[gene]    = 1
         gnreftx[gene]  = RefTx[ RefTx( txinfo, trandon, tranacc, txlen ) ]
      end
   end
   # now make RefSet and add genes.
   for gene in keys(gninfo)
      geneset[gene] = RefGene( gninfo[gene],  
                               gndon[gene],   gnacc[gene],
                               gntxst[gene],  gntxen[gene], 
                               gnexons[gene], 
                               gnlen[gene] /  gncnt[gene], 
                               gnreftx[gene] )
   end

   return geneset
end

function fetch_meta( var::String, meta; off=1 )
   val = "_" # empty
   for i in off:2:length(meta)
      if meta[i] == var && i <= length(meta)
         return string(meta[i+1]),i+1
      end
   end
   val,0
end

function load_gtf( fh; txbool=true, suppress=false )

   # Temporary variables   
   gninfo   = Dict{GeneName,GeneInfo}()
   gndon    = Dict{GeneName,CoordTuple}()
   gnacc    = Dict{GeneName,CoordTuple}()
   gntxst   = Dict{GeneName,CoordTuple}()
   gntxen   = Dict{GeneName,CoordTuple}()
   gnlens   = Dict{GeneName,CoordTuple}()
   gnexons  = Dict{GeneName,CoordTree}()
   gnreftx  = Dict{GeneName,Vector{RefTx}}()
   gnlen    = Dict{GeneName,Float64}()
   gncnt    = Dict{GeneName,Int}()

   # RefSet variables
   geneset  = Dict{GeneName,RefGene}()

   used_txn = Dict{GeneName,Bool}()
   warning  = false

   txnum = 75000
   sizehint!(geneset, txnum >> 1)

   curtran = ""
   curgene = ""
   curchrom = ""
   curstran = '+'

   trandon = Vector{CoordInt}()
   tranacc = Vector{CoordInt}()
   txlen   = 0

   function private_add_transcript!( curtran, curgene, curchrom, curstran, trandon, tranacc, txlen )
      sort!(trandon)
      sort!(tranacc)

      txS = minimum(tranacc)
      txE = maximum(trandon)

      txinfo = TxInfo(curtran, txS, txE, UInt16(length(tranacc)))

      don = tuple(trandon[1:(end-1)]...)
      acc = tuple(tranacc[2:end]...)

      if haskey(gninfo, curgene)
         curchrom == gninfo[curgene].name || return # can have only one chrom
         gndon[curgene]  = unique_tuple(gndon[curgene], don)
         gnacc[curgene]  = unique_tuple(gnacc[curgene], acc)
         gntxst[curgene] = unique_tuple(gntxst[curgene], tuple(txS))
         gntxen[curgene] = unique_tuple(gntxen[curgene], tuple(txE))
         gnlen[curgene] += txlen
         gncnt[curgene] += 1
         push!( gnreftx[curgene], RefTx( txinfo, tuple(trandon...), tuple(tranacc...), txlen ) )
      else
         gndon[curgene]  = don
         gnacc[curgene]  = acc
         gninfo[curgene] = GeneInfo(curgene, string(curchrom), curstran)
         gntxst[curgene] = tuple(txS)
         gntxen[curgene] = tuple(txE)
         gnlen[curgene]  = txlen
         gncnt[curgene]  = 1
         gnreftx[curgene]  = RefTx[ RefTx( txinfo, tuple(trandon...), tuple(tranacc...), txlen ) ]
      end
   end

   for l in eachline(fh)

      (l[1] == '#') && continue # ignore comment lines

      (chrom, src, entrytype,
       st, en, _, 
       strand, _,
       meta) = split(chomp(l), '\t')

      (entrytype == "exon") || continue

      metaspl = split(meta, [' ',';','"'], keep=false)

      geneid,pos  = fetch_meta( "gene_id", metaspl )
      tranid,pos  = fetch_meta( "transcript_id", metaspl )
      genesym     = fetch_meta( "gene_name", metaspl )
      support,val = fetch_meta( "transcript_support_level", metaspl )

      # if we observe low transcript support levels, we will warn the user
      if (support != "_" && support != "1" && support != "2" &&  support != "NA")

         if suppress
            continue
         elseif !warning
            println(STDERR, "")
            warn("Using low quality Transcript Support Levels (TSL 3+) in your GTF file is not recommended!\nFor more information on TSL, see: http://www.ensembl.org/Help/Glossary?id=492\n\nIf you would like Whippet to ignore these when building its index, use `--suppress-low-tsl` option!\n\n")
         
            warning=true
         end

      elseif haskey(used_txn, tranid)
         error("GTF file is not in valid GTF2.2 format!\n\nAnnotation entries for 'transcript_id' $tranid has already been fully processed and closed.\nHint: All GTF lines with the same 'transcript_id' must be adjacent in the GTF file and referring to the same transcript and gene!")
      elseif tranid == geneid && tranid != curtran
         warn("'transcript_id' should not equal 'gene_id' at $tranid == $geneid")
      end

      if tranid != curtran
         # add transcript and clean up
         curtran != "" && private_add_transcript!( curtran, curgene, curchrom, curstran, trandon, tranacc, txlen )

         used_txn[curtran] = true

         curtran = tranid
         curgene = geneid
         curchrom = chrom
         curstran = strand[1]

         empty!(trandon)
         empty!(tranacc)
         txlen = 0

      end 

      sval = parse_coordint(st)
      eval = parse_coordint(en)

      push!(tranacc, sval) 
      push!(trandon, eval)
      txlen += eval - sval + 1

      insval = Interval{CoordInt}(sval, eval)
      if haskey(gnexons, curgene)
         # make sure we are adding a unique value
         if !haskey(gnexons[curgene], (sval, eval))
            push!(gnexons[curgene], insval)
         end
      else
         gnexons[curgene] = CoordTree()
         push!(gnexons[curgene], insval)
      end 
   end

   private_add_transcript!( curtran, curgene, curchrom, curstran, trandon, tranacc, txlen )
   
   for gene in keys(gninfo)
      geneset[gene] = RefGene( gninfo[gene],
                               gndon[gene],   gnacc[gene],
                               gntxst[gene],  gntxen[gene],
                               gnexons[gene],
                               gnlen[gene] /  gncnt[gene],
                               gnreftx[gene] )
   end

   return geneset
end

