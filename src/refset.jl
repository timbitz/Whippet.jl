
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
   noveldon::CoordTuple
   novelacc::CoordTuple
   exonexpr::Float64
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
unique_tuple( tup1::Tuple{Vararg{T}}, tup2::Tuple{} ) where T = unique_tuple(tup1, tup1)
unique_tuple( tup1::Tuple{}, tup2::Tuple{Vararg{T}} ) where T = unique_tuple(tup2, tup2)

function unique_tuple( tup1::Tuple{Vararg{T}}, tup2::Tuple{Vararg{T}}) where T
   uniq = SortedSet{T,Base.Order.ForwardOrdering}()
   for i in tup1, j in tup2
     push!(uniq, i)
     push!(uniq, j)
   end
   tuple(uniq...)
end

function minimum_threshold!( dict::Dict{K,V}, limit::V, range::UnitRange ) where {K, V <: Number}
   for k in collect(keys(dict))
      if dict[k] < limit || !(k in range)
         delete!( dict, k )
      end
   end
end

function load_refflat( fh; txbool=true )
   error("ERROR: Annotation files in refflat format have been deprecated as of Whippet v0.11, please use --gtf!")
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

function load_gtf( fh; txbool=true, suppress=false, usebam=false, bamreader=Nullable{BAM.Reader}(), bamreads=2, bamoneknown=false )

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
   warning_num = 0

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

      metaspl = split(meta, [' ',';','"'], keepempty=false)

      geneid,pos  = fetch_meta( "gene_id", metaspl )
      tranid,pos  = fetch_meta( "transcript_id", metaspl )
      genesym     = fetch_meta( "gene_name", metaspl )
      support,val = fetch_meta( "transcript_support_level", metaspl )

      # if we observe low transcript support levels, we will warn the user
      if (support != "_" && support != "1" && support != "2" &&  support != "NA")

         if suppress
            continue
         elseif !warning
            println(stderr, "")
            @warn("Using low quality Transcript Support Levels (TSL 3+) in your GTF file is not recommended!\nFor more information on TSL, see: http://www.ensembl.org/Help/Glossary?id=492\n\nIf you would like Whippet to ignore these when building its index, use `--suppress-low-tsl` option!\n\n")

            warning=true
         end

      elseif haskey(used_txn, tranid)
         error("ERROR: GTF file is not in valid GTF2.2 format!\n\nERROR: Annotation entries for 'transcript_id' $tranid has already been fully processed and closed.\nHINT: All GTF lines with the same 'transcript_id' must be adjacent in the GTF file and referring to the same transcript and gene!")
      elseif tranid == geneid && tranid != curtran
         if warning_num < 25
            @warn("Generally 'transcript_id' should not equal 'gene_id' but does at $tranid == $geneid;")
         elseif warning_num == 25
            @warn("... NOTE: 'transcript_id' == 'gene_id' will work OK for single isoform genes, but will not produce expected behavior for multi-isoform genes!\n... similar warnings will be suppressed; disregard if 'transcript_id' == 'gene_id' is intentional\n")
         end
         warning_num += 1
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

      insval = IntervalTrees.Interval{CoordInt}(sval, eval)
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

   if usebam
      bamnames = map(x->values(x)[1], findall(BAM.header(get(bamreader)), "SQ"))
   end

   novel_ss = 0
   annot_ss = 0
   # iterate through genes, if usebam, add novel splice sites and calculate relative exonexpr
   # set RefGene at the end of each iteration
   for gene in keys(gninfo)

      novelacc = Dict{CoordInt,Int}()
      noveldon = Dict{CoordInt,Int}()
      exonexpr = 0.0
      meanleng = gnlen[gene] / gncnt[gene]
      seqname  = gninfo[gene].name
      strand   = gninfo[gene].strand
      bamwasused = false

      annot_ss += length(gndon[gene]) + length(gnacc[gene])

      if usebam && seqname in bamnames
         bamwasused = true
         range   = Int64(minimum(gntxst[gene])):Int64(maximum(gntxen[gene]))
         exoncount = process_records!( get(bamreader), seqname, range, strand,
                                       gnexons[gene], union(gnacc[gene], gndon[gene]), bamoneknown,
                                       novelacc, noveldon )
         exonexpr  = exoncount / (meanleng - 100)

         # clean up cryptic splice-sites
         minimum_threshold!( novelacc, bamreads, range )
         minimum_threshold!( noveldon, bamreads, range )
         novelacctup = map(CoordInt, Tuple(collect(keys(novelacc))))
         noveldontup = map(CoordInt, Tuple(collect(keys(noveldon))))

         # add splice-sites from bam
         annoacc     = gnacc[gene]
         annodon     = gndon[gene]
         gnacc[gene] = unique_tuple( gnacc[gene], novelacctup )
         gndon[gene] = unique_tuple( gndon[gene], noveldontup )
         novel_ss += length(gnacc[gene]) + length(gndon[gene]) - (length(annoacc) + length(annodon))
      end

      # clean up any txst or txen that match a splice site
      gntxst[gene] = Tuple( setdiff( gntxst[gene], gnacc[gene] ) )
      gntxen[gene] = Tuple( setdiff( gntxen[gene], gndon[gene] ) )

      geneset[gene] = RefGene( gninfo[gene],
                               gndon[gene],
                               gnacc[gene],
                               gntxst[gene],
                               gntxen[gene],
                               gnexons[gene],
                               bamwasused ? Tuple( setdiff(noveldontup, annodon) ) : Tuple{}(),
                               bamwasused ? Tuple( setdiff(novelacctup, annoacc) ) : Tuple{}(),
                               exonexpr,
                               meanleng,
                               gnreftx[gene] )
   end

   if usebam
      println(stderr, "Found $novel_ss new splice-sites from BAM file..")
   end
   println(stderr, "Loaded $annot_ss annotated splice-sites from GTF file..")

   return geneset
end
