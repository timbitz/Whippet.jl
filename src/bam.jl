
const AlignBlocks = Vector{Tuple{CoordInt,CoordInt}}
const AlignNodes  = Vector{GenomicFeatures.Interval{SGNodeIsExon}}

const EMPTY_PATH   = Vector{SGNodeIsExon}()
const SPLICE_BONUS = 10
const MIN_OVERLAP  = 6

# reuse temporary data
struct AlignData
   blocks::AlignBlocks
   nodes::AlignNodes
end

function AlignData()
   blocks   = AlignBlocks()
   nodes    = AlignNodes()
   sizehint!(blocks, 100)
   sizehint!(nodes,  100)
   AlignData(blocks, nodes)
end

function Base.empty!( data::AlignData )
   empty!( data.blocks )
   empty!( data.nodes )
   data
end

tempname( rec::SAM.Record ) = SAM.tempname(rec)
tempname( rec::BAM.Record ) = BAM.tempname(rec)
refname( rec::SAM.Record ) = SAM.refname(rec)
refname( rec::BAM.Record ) = BAM.refname(rec)
nextrefname( rec::SAM.Record ) = SAM.nextrefname(rec)
nextrefname( rec::BAM.Record ) = BAM.nextrefname(rec)
flag( rec::SAM.Record ) = SAM.flag(rec)
flag( rec::BAM.Record ) = BAM.flag(rec)
leftposition( rec::SAM.Record ) = SAM.leftposition(rec)
leftposition( rec::BAM.Record ) = BAM.leftposition(rec)
rightposition( rec::SAM.Record ) = SAM.rightposition(rec)
rightposition( rec::BAM.Record ) = BAM.rightposition(rec)
sequence( rec::SAM.Record ) = SAM.sequence(rec)
sequence( rec::BAM.Record ) = BAM.sequence(rec)
hasalignment( rec::SAM.Record ) = SAM.hasalignment(rec)
hasalignment( rec::BAM.Record ) = BAM.hasalignment(rec)
hasposition( rec::SAM.Record ) = SAM.hasposition(rec)
hasposition( rec::BAM.Record ) = BAM.hasposition(rec)

ismapped( rec::SAM.Record ) = SAM.ismapped(rec)
ismapped( rec::BAM.Record ) = BAM.ismapped(rec)
isnextmapped( rec::SAM.Record ) = SAM.isnextmapped(rec)
isnextmapped( rec::BAM.Record ) = BAM.isnextmapped(rec)

cigar( rec::SAM.Record ) = SAM.cigar(rec)
cigar( rec::BAM.Record ) = BAM.cigar(rec)
cigar_rle( rec::BAM.Record ) = BAM.cigar_rle(rec)

function cigar_rle( rec::SAM.Record )
   opstrings = split( cigar(rec), r"(?<=[IDMNXSH=])(?=\d)" )
   op_arr  = Vector{BioAlignments.Operation}(undef, length(opstrings))
   len_arr = Vector{Int64}(undef, length(opstrings))
   for i in 1:length(opstrings)
      cur_op = split( opstrings[i], r"(?<=\d)(?=[IDMNXSH=])" )
      len_arr[i] = parse(Int, cur_op[1])
      op_arr[i] = convert(BioAlignments.Operation, cur_op[2][1])
   end
   (op_arr, len_arr)
end

strandchar( rec::R ) where R <: Union{SAM.Record, BAM.Record} = isstrandfwd(rec) ? '+' : '-'
revstrandchar( rec::R ) where R <: Union{SAM.Record, BAM.Record} = isstrandfwd(rec) ? '-' : '+'
isskipchar( c::Char ) = c == 'N'
introncount( rec::R ) where R <: Union{SAM.Record, BAM.Record} = count(isskipchar, cigar(rec))
isspliced( rec::R ) where R <: Union{SAM.Record, BAM.Record} = introncount(rec) > 0
isstrandfwd( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x010 == 0
ispaired( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x1 != 0
isleftmate( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x40 != 0
isrightmate( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x80 != 0

const FLOAT_SHIFT  = 100_000_000
const MAX_POSITION = 9_999_999

function encode_aberrant( node, pos )
   @assert pos <= MAX_POSITION # maximum position offset in node
   NodeFloat(node + pos / FLOAT_SHIFT)
end

function decode_aberrant( enc::NodeFloat )
   root = NodeInt( floor(enc) )
   node = trunc(enc, digits=1)
   pos  = Int(round(Float64(enc - node) * FLOAT_SHIFT))
   root, node, pos
end

function decode_aberrant( enc::I ) where I <: Integer
   NodeInt(enc), NodeFloat(enc), 0
end

isexonic( node::NodeInt )   = false
function isexonic( node::NodeFloat )
   root, node, pos = decode_aberrant( node )
   (root == node) ? true : false
end

function flipaberrant!( v::Vector{SGNodeIsExon} )
   for i in 1:length(v)
      n = v[i]
      root,node,pos = decode_aberrant( n.node )
      if isaberrant(n.node) && pos >= 1
         v[i] = SGNodeIsExon(n.gene, n.node, !n.meta)
      end
   end
end

function fill_alignment_blocks!( blocks::AlignBlocks, 
                                 rec::R ) where R <: Union{SAM.Record, BAM.Record}
   op,len = cigar_rle( rec )
   left   = leftposition(rec)
   right  = rightposition(rec)
   refpos = left
   empty!(blocks)

   for i in 1:length(op)
      if ismatchop(op[i])  
         refpos += len[i]
      elseif isdeleteop(op[i])
         donor = CoordInt(refpos - 1)
         refpos += len[i]
         accep = CoordInt(refpos)
         op[i] == OP_SKIP || continue

         push!( blocks, (left, donor))
         left = accep
      end
   end
   push!( blocks, (left, right))
   blocks
end

#= deprecated
function extract_alignment_blocks( rec::R ) where R <: Union{SAM.Record, BAM.Record}
   blocks = AlignBlocks()
   fill_alignment_blocks!( blocks, rec )
   blocks
end


function extract_junctions( rec::R ) where R <: Union{SAM.Record, BAM.Record}
   op,len = cigar_rle( rec )
   refpos = leftposition( rec )

   introns = Vector{Tuple{CoordInt,CoordInt}}()

   for i in 1:length(op)
      if ismatchop(op[i])  
         refpos += len[i]
      elseif isdeleteop(op[i])
         donor = CoordInt(refpos - 1)
         refpos += len[i]
         accep = CoordInt(refpos)
         op[i] == OP_SKIP || continue
 
         push!(introns, (donor, accep))
      end
   end
   introns
end=#

function unstranded_overlap( i_a::I, i_b::I,
                             j_a::J, j_b::J ) where {I <: Integer, J <: Integer}
   max(0, min(i_b, j_b) - max(i_a, j_a) + 1)
end

# return key that has the first max value
function max_value_key( d::Dict{K,V} ) where {K,V}
   bestval = first(values(d))
   bestkey = first(keys(d))
   for (k,v) in d
      if v > bestval
         bestkey = k
         bestval = v
      end
   end
   bestkey
end

function overlapping_nodes!( ic::IntervalCollection{SGNodeIsExon}, 
                             data::AlignData,
                             overlaps::Dict{NodeInt, Int},
                             rec::R;
                             flipstrand=false ) where R <: Union{SAM.Record, BAM.Record}
   
   refn   = refname(rec)
   stranded = length(data.blocks) > 1 ? true : false

   for i in 1:length(data.blocks)
      l,r = data.blocks[i]
      
      for n in eachoverlap( ic, GenomicFeatures.Interval(refn, l, r))
         if !stranded || rec["XS"] == convert(Char, n.strand)
            over = unstranded_overlap( n.first, n.last, l, r)
            over = n.metadata.meta ? over : Int(floor(over/2))
            over += (i > 1 && n.last == r) ? SPLICE_BONUS : 0
            over += (i <= length(data.blocks) && n.first == l) ? SPLICE_BONUS : 0
            increment!( overlaps, n.metadata.gene, over)
            push!( data.nodes, n)
         end
      end
   end
end

function best_path( data::AlignData, bestgene::GeneInt, min_over::Int=MIN_OVERLAP )

   path     = Vector{SGNodeIsExon}()
   inde     = Vector{Int}()
   sizehint!(path, length(data.nodes))

   # pull all the overlapping nodes
   for i in 1:length(data.nodes)
      if data.nodes[i].metadata.gene == bestgene
         push!( path, data.nodes[i].metadata )
         push!( inde, i )
      end
   end

   # check tails to ensure minimum nucleotide overlap
   if length(path) > 1
      # check overlap of first node
      leng = data.nodes[inde[1]].last - data.nodes[inde[1]].first + 1
      over = unstranded_overlap( data.nodes[inde[1]].first, 
                                 data.nodes[inde[1]].last, 
                                 data.blocks[1][1], 
                                 data.blocks[1][2] )
      if over < leng && over < min_over
         popfirst!(path)
      end

      # check last node
      leng = data.nodes[inde[length(inde)]].last - data.nodes[inde[length(inde)]].first + 1
      over = unstranded_overlap( data.nodes[inde[length(inde)]].first, 
                                 data.nodes[inde[length(inde)]].last, 
                                 data.blocks[length(data.blocks)][1], 
                                 data.blocks[length(data.blocks)][2] )
      if over < leng && over < min_over
         pop!(path)
      end
   end

   path
end


function all_introns_valid( data::AlignData, gene::GeneInt )
   j = 1
   for i in 1:length(data.blocks)
      if i > 1
         for k in j:length(data.nodes)
            j = k
            data.nodes[k].metadata.gene == gene || continue
            if first(data.nodes[k]) == data.blocks[i][1]
               break
            elseif first(data.nodes[k]) > data.blocks[i][1]
               return false
            end
         end
      end

      if i < length(data.blocks)
         for k in j:length(data.nodes)
            j = k
            data.nodes[k].metadata.gene == gene || continue
            if last(data.nodes[k]) == data.blocks[i][2]
               j += 1
               break
            elseif last(data.nodes[k]) > data.blocks[i][2]
               return false  
            end
         end
      end
   end
   true
end

function aberrant_path( data::AlignData, 
                        gene::GeneInt, 
                        strand::Bool=true, 
                        verbose=false )

   path     = Vector{SGNodeIsExon}()
   sizehint!(path, length(data.nodes))

   j = 1
   for i in 1:length(data.blocks)

      if i > 1
         for k in j:length(data.nodes)
            j = k
            data.nodes[k].metadata.gene == gene || continue

            if first(data.nodes[k]) == data.blocks[i][1]
               push!( path, data.nodes[k].metadata )
               break
            elseif first(data.nodes[k]) < data.blocks[i][1] <= last(data.nodes[k])
               pos  = strand ? data.blocks[i][1] - first(data.nodes[k]) : last(data.nodes[k]) - data.blocks[i][1]
               aber = encode_aberrant( data.nodes[k].metadata.node, pos + 1)
               span = SGNodeIsExon( data.nodes[k].metadata.gene, aber :: NodeNum, false)
               push!( path, span )
               break
            end
         end 
      end

      if i < length(data.blocks)
         for k in j:length(data.nodes)
            j = k
            data.nodes[k].metadata.gene == gene || continue

            if last(data.nodes[k]) == data.blocks[i][2]
               data.nodes[k].metadata in path || push!( path, data.nodes[k].metadata )
               j += 1
               break
            elseif first(data.nodes[k]) <= data.blocks[i][2] < last(data.nodes[k])
               pos  = strand ? data.blocks[i][2] - first(data.nodes[k]) : last(data.nodes[k]) - data.blocks[i][2]
               aber = encode_aberrant( data.nodes[k].metadata.node, pos + 1)
               span = SGNodeIsExon( data.nodes[k].metadata.gene, aber :: NodeNum, true)
               push!( path, span )
               break
            elseif data.blocks[i][1] <= first(data.nodes[k]) &&
                   last(data.nodes[k])  <  data.blocks[i][2]
               push!( path, data.nodes[k].metadata )
            end
         end
      end
   end
   verbose && println(stderr, path)
   # add canonical node to remove dangling aberrant edges
   if length(path) > 1 && isaberrant(path[1].node)
      root,node,pos = decode_aberrant(path[1].node)
      if root == node
         pushfirst!(path, SGNodeIsExon(path[1].gene, root, true))
      else
         pushfirst!(path, SGNodeIsExon(path[1].gene, node, false))
      end
   end
   if length(path) > 1 && isaberrant(path[length(path)].node)
     root,node,pos = decode_aberrant(path[length(path)].node)
      if root == node
         push!(path, SGNodeIsExon(path[length(path)].gene, root, true))
      else
         push!(path, SGNodeIsExon(path[length(path)].gene, node, false))
      end
   end
   path
end

function choose_path( data::AlignData, 
                      gene::GeneInt, 
                      strand::Bool=true, 
                      allow_aberrant::Bool=true )

   if allow_aberrant && !all_introns_valid( data, gene )
      return aberrant_path( data, gene, strand ), false
   else
      return best_path( data, gene ), true
   end
end

function align_bam( ic::IntervalCollection{SGNodeIsExon},
                    info::Vector{GeneInfo},
                    data::AlignData,
                    overlaps::Dict{NodeInt, Int},
                    rec::R,
                    allow_aberrant::Bool=true ) where R <: Union{SAM.Record, BAM.Record}
   empty!(data) 
   empty!(overlaps)
   fill_alignment_blocks!( data.blocks, rec )
   overlapping_nodes!(ic, data, overlaps, rec)
   if length(overlaps) > 0 
      gene = max_value_key(overlaps)
      path, aber = choose_path( data, gene, info[gene].strand, allow_aberrant )

      if length(path) > 0
         return path, aber
      end
   end
   EMPTY_PATH, false
end

function align_bam( ic::IntervalCollection{SGNodeIsExon},
                    info::Vector{GeneInfo},
                    leftdata::AlignData,
                    rightdata::AlignData,
                    overlaps::Dict{NodeInt, Int},
                    leftmate::R,
                    rightmate::R,
                    allow_aberrant::Bool=true ) where R <: Union{SAM.Record, BAM.Record}
   empty!(leftdata) 
   empty!(rightdata)
   empty!(overlaps)
   fill_alignment_blocks!( leftdata.blocks, leftmate )
   fill_alignment_blocks!( rightdata.blocks, rightmate )
   overlapping_nodes!(ic, leftdata, overlaps, leftmate)
   overlapping_nodes!(ic, rightdata, overlaps, rightmate, flipstrand=true)
   if length(overlaps) > 0 
      gene      = max_value_key(overlaps)
      leftpath, laber  = choose_path( leftdata, gene, info[gene].strand, allow_aberrant )
      rightpath,raber  = choose_path( rightdata, gene, info[gene].strand, allow_aberrant )
      #println(stderr, "left strand: $(isstrandfwd(leftmate))")
      #println(stderr, leftpath)
      #println(stderr, all_introns_valid( leftdata, gene ))
      #println(stderr, "right strand: $(isstrandfwd(rightmate))")
      #println(stderr, rightpath)
      #println(stderr, all_introns_valid( rightdata, gene ))
      #println(stderr, "gene: $gene")
      #println(stderr, lib.info[gene])
      if length(leftpath) > 0 && 
         length(rightpath) > 0 
         return leftpath, laber, 
                rightpath, raber
      end
   end
   EMPTY_PATH, false, 
   EMPTY_PATH, false
end

function make_bamparser( bamfile )
   isfile(bamfile) || error("ERROR: --bam parameter used, but cannot find .bam file at $bamfile !")
   isfile(bamfile * ".bai") || error("ERROR: no .bai index found for $bamfile! Cannot set-up random access to $bamfile !")
   println(stderr, "Loading BAM file for random-access: $bamfile")
   bamreadr = open(BAM.Reader, bamfile, index=bamfile * ".bai")
   bamreadr
end

function ispaired_bamfile( bamfile, max_reads=1000 )
   parser = open(BAM.Reader, bamfile)
   i = 0
   while !eof(parser) && i < max_reads
      ispaired(first(parser)) && (close(parser); return true)
      i += 1
   end
   close(parser)
   return false
end

function has_collated_header( parser, so_value="unsorted", go_value="query" )
   fhead = findall(SAM.header(parser), "HD")
   if length(fhead) == 1
      header_order = fhead[1]
      header_order["SO"] == so_value || return false
      header_order["GO"] == go_value || return false
      return true
   else
      return false
   end
end

function check_samtools()
   hndl = open(`samtools --version`, "r")
   readline(hndl)[1:8] == "samtools" || (close(hndl) && return false)
   close(hndl)
   hndl = open(pipeline(`samtools --help`, `grep collate`), "r")
   collate = readlines(hndl)
   close(hndl)
   if length(collate) == 1
      return true
   else
      return false
   end
end

function make_single_parser( bamfile )
   isfile(bamfile) || error("ERROR: --bam parameter used, but cannot find .bam file at $bamfile !")
   println(stderr, "Loading BAM file single-end reads...")
   parser = open( BAM.Reader, bamfile )
   parser
end

function make_collated_parser( filename; tmpdir=ENV["TMPDIR"], flag="-f", force_collated=true )
   check_samtools() || error("Samtools is not properly installed, but is a dependency of BAM input!")
   prefix = rstrip(tmpdir, ['/']) * "/" * splitext(basename(filename))[1]
   head = SAM.Reader(open(`samtools view -H $filename`, "r"))
   iscollated = has_collated_header(head)

   if !iscollated || force_collated
      println(stderr, "Loading BAM file and collating with samtools...")
      println(stderr, "TMPDIR: " * tmpdir)
      @time parser = SAM.Reader(open(`samtools collate $filename $prefix --output-fmt SAM -O $flag`, "r"))
      has_collated_header(parser) || error("samtools collate failed to produce a proper input file!")
      iscollated = true
   else
      if iscollated
         println(stderr, "Loading header verified collated BAM file...")
      else
         println(stderr, "Warning: BAM file is not collated! ")
         println(stderr, "Loading BAM in single-end primary-only mode")
      end
      @time parser = open(BAM.Reader, filename)
   end
   parser #, iscollated
end



function read_collated_bam( parser, 
                            lib,
                            quant,
                            mod;
                            stranded=false,
                            fr_oriented=true,
                            allow_aberrant=true )
   concordant    = 0
   nonconcordant = 0
   singletons    = 0
   numintrons    = 0

   mean_readlen  = 0.0
  
   leftmate  = eltype(parser)()
   rightmate = eltype(parser)()

   maxdist = 5000000

   i = 1

   leftdata  = AlignData()
   rightdata = AlignData()
   overlaps  = Dict{NodeInt, Int}()
   sizehint!(overlaps, 16)

   while !eof( parser )

      read!( parser, leftmate )
      read!( parser, rightmate )
      
      while tempname(leftmate) != tempname(rightmate)
         # Process leftmate as singleton
         #align and count

         if hasalignment(leftmate)
            leftpath, leftvalid = align_bam( lib.coords, 
                                             lib.info, 
                                             leftdata, 
                                             overlaps, 
                                             leftmate, 
                                             allow_aberrant )
            if length(leftpath) > 0 

               leftseq = sequence(leftmate) 
               if !isstrandfwd(leftmate)
                  BioSequences.reverse_complement!(leftseq)
               end
               biasval = count!( mod, leftseq )

               if !lib.info[leftpath[1].gene].strand
                  reverse!(leftpath)
                  !leftvalid && flipaberrant!(leftpath)
               end

               count!( quant, leftpath, biasval )
               @fastmath mean_readlen += (length(leftseq) - mean_readlen) / i
               i += 1
            end
         end

         singletons += 1

         # Now shift
         leftmate,rightmate = rightmate,leftmate
         empty!(rightmate)
         eof(parser) && break
         read!(parser, rightmate)
         #println("shifting singleton at $i")
      end
      eof(parser) && break

      hasalignment(leftmate) && hasalignment(rightmate) || continue

      dist = leftposition(rightmate) - leftposition(leftmate)
      if !isleftmate(leftmate) || 
         !isrightmate(rightmate) ||
         refname(leftmate) != refname(rightmate) ||
         abs(dist) > maxdist ||
         (isstrandfwd(leftmate) && dist < 0) ||
         (!isstrandfwd(leftmate) && dist > 0)


         # Non-concordant reads
         nonconcordant += 2

      else

         #println(stderr, "dist is: $dist > 0? && left: $(strandchar(leftmate))")
         # Process left/right first_mates
         leftpath, leftvalid, rightpath, rightvalid  = align_bam( lib.coords,
                                                                  lib.info, 
                                                                  leftdata, 
                                                                  rightdata,
                                                                  overlaps,
                                                                  leftmate,
                                                                  rightmate,
                                                                  allow_aberrant )

         #println(stderr, "left: $leftvalid, right: $rightvalid")

         if length(leftpath) > 0  &&
            length(rightpath) > 0 

            leftseq  = sequence(leftmate)
            rightseq = sequence(rightmate)

            if isstrandfwd(leftmate)
               biasval = count!( mod, leftseq, rightseq )
            else
               BioSequences.reverse_complement!(leftseq)
               BioSequences.reverse_complement!(rightseq)
               biasval = count!( mod, leftseq, rightseq )
            end

            if !lib.info[leftpath[1].gene].strand
               reverse!(leftpath)
               reverse!(rightpath)
               !leftvalid && flipaberrant!(leftpath)
               !rightvalid && flipaberrant!(rightpath)
            end

            if (dist < 0 && lib.info[leftpath[1].gene].strand) ||
               (dist > 0 && !lib.info[leftpath[1].gene].strand)
               leftpath, rightpath = rightpath, leftpath
            end

            #if length(leftpath) > 1 || length(rightpath) > 1
                  #println(stderr, "$(isstrandfwd(leftmate)) != $(lib.info[leftpath[1].gene].strand)")
                  #println(stderr, leftpath)
                  #println(stderr, rightpath)
            #end

            count!( quant, leftpath, rightpath, biasval )

            @fastmath mean_readlen += (length(leftseq) - mean_readlen) / i
            i += 1
            @fastmath mean_readlen += (length(rightseq) - mean_readlen) / i
            i += 1

            #println(stderr, leftmate)
            #println(stderr, isstrandfwd(leftmate))
            #println(stderr, isleftmate(leftmate))
            #println(stderr, leftpath)
            #println(stderr, "")
            #println(stderr, rightmate)
            #println(stderr, isstrandfwd(rightmate))
            #println(stderr, isrightmate(rightmate))
            #println(stderr, rightpath)
            #error()
         end

         concordant += 2
      end

      if i % 10000000 == 0
         GC.gc()
      end

      empty!(leftmate)
      empty!(rightmate)
   end
   concordant, nonconcordant, singletons, mean_readlen
end





function process_records!( reader::BAM.Reader, 
                           seqname::String, 
                           range::UnitRange{Int64}, 
                           strand::Bool,
                           exons::CoordTree, 
                           known, 
                           oneknown::Bool,
                           novelacc::Dict{K,V}, 
                           noveldon::Dict{K,V} ) where {K,V}
   exoncount = 0
   try
      for rec in GenomicFeatures.eachoverlap( reader, seqname, range )
         # add to exon count if overlaps exon
         if hasintersection( exons, leftposition(rec) ) ||
            hasintersection( exons, rightposition(rec) )
            exoncount += 1
         end
         # if is spliced process splice sites
         if isspliced(rec) && strand == isstrandfwd(rec)
            known = process_spliced_record!( novelacc, noveldon, rec, known, oneknown )
         end
      end
   catch e
      error(e)
   end
   exoncount
end

function process_spliced_record!( novelacc::Dict{K,V}, 
                                  noveldon::Dict{K,V}, 
                                  rec::BAM.Record,
                                  known, 
                                  oneknown::Bool ) where {K,V}
   op,len = BAM.cigar_rle( rec )
   refpos = BAM.position( rec )
   for i in 1:length(op)
      if ismatchop(op[i])  
         refpos += len[i]
      elseif isdeleteop(op[i])
         donor = CoordInt(refpos - 1)
         refpos += len[i]
         accep = CoordInt(refpos)
         op[i] == OP_SKIP || continue
         if oneknown && !(donor in known) && !(accep in known)
            continue
         else
            increment!( noveldon, donor, 1 )
            increment!( novelacc, accep, 1 )
            known = union(known, (donor, accep))
         end
      end
   end
   known
end