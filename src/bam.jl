
const AlignBlocks = Vector{Tuple{CoordInt,CoordInt}}
const AlignNodes  = Vector{GenomicFeatures.Interval{SGNodeIsExon}}

const EMPTY_PATH   = Vector{SGNodeIsExon}()

# reuse temporary data
struct AlignData
   blocks::AlignBlocks
   nodes::AlignNodes
   overlaps::Dict{NodeInt, Int}
end

function AlignData()
   blocks   = AlignBlocks()
   nodes    = AlignNodes()
   overlaps = Dict{NodeInt,Int}()
   sizehint!(blocks, 100)
   sizehint!(nodes,  100)
   sizehint!(overlaps, 25)
   AlignData(blocks, nodes, overlaps)
end

function Base.empty!( data::AlignData )
   empty!( data.blocks )
   empty!( data.overlaps )
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
isskipchar( c::Char ) = c == 'N'
introncount( rec::R ) where R <: Union{SAM.Record, BAM.Record} = count(isskipchar, cigar(rec))
isspliced( rec::R ) where R <: Union{SAM.Record, BAM.Record} = introncount(rec) > 0
isstrandfwd( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x010 == 0
ispaired( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x1 != 0
isleftmate( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x40 != 0
isrightmate( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x80 != 0


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
   max(0, min(i_b, j_b) - max(i_a, j_a))
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

function overlapping_nodes( ic::IntervalCollection{SGNodeIsExon}, 
                            data::AlignData,
                            rec::R ) where R <: Union{SAM.Record, BAM.Record}
   
   refn   = refname(rec)

   for i in 1:length(data.blocks)
      l,r = data.blocks[i]
      
      for n in eachoverlap( ic, GenomicFeatures.Interval(refn, l, r))
         if strandchar(rec) == convert(Char, n.strand)
            overlap = unstranded_overlap(n.first, n.last, l, r)
            increment!( data.overlaps, n.metadata.gene, overlap )
            push!( data.nodes, n)
         end
      end
   end
   length(data.overlaps) > 0 || return EMPTY_PATH, zero(NodeInt)

   bestgene = max_value_key(data.overlaps)
   path     = Vector{SGNodeIsExon}()
   sizehint!(path, length(data.nodes))

   for n in data.nodes
      if n.metadata.gene == bestgene
         push!( path, n.metadata )
      end
   end

   path, bestgene
end

function check_introns( data )
   j = 1
   badsplice = Vector{eltype(eltype(data.blocks))}()
   for i in 1:length(data.blocks)
      println(Int(data.blocks[i][1]))
      println(Int(data.blocks[i][2]))
      println("$i\n")
      if i > 1
         for k in j:length(data.nodes)
            println(k)
            println(data.nodes[k])
            j = k
            if first(data.nodes[k]) == data.blocks[i][1]
               println(first(data.nodes[k]))
               println(Int(data.blocks[i][1]))
               break
            elseif first(data.nodes[k]) > data.blocks[i][1]
               push!(badsplice, data.blocks[i][1])
               break
            end
         end
      end
      println()
      if i < length(data.blocks)
         for k in j:length(data.nodes)
            println(k)
            println(data.nodes[k])
            j = k
            if last(data.nodes[k]) == data.blocks[i][2]
               println(last(data.nodes[k]))
               println(Int(data.blocks[i][2]))
               j += 1
               break
            elseif last(data.nodes[k]) > data.blocks[i][2]
               push!(badsplice, data.blocks[i][2])
               break
            end
         end
      end
      println()
   end
   badsplice
end

function all_introns_valid( data::AlignData, gene::NodeInt )
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


function align_bam( ic::IntervalCollection{SGNodeIsExon}, 
                    data::AlignData,
                    rec::R ) where R <: Union{SAM.Record, BAM.Record}
   empty!(data) 
   fill_alignment_blocks!( data.blocks, rec )
   path,gene = overlapping_nodes(ic, data, rec)
   if length(path) > 0 && all_introns_valid( data, gene )
      return path, true
   end
   EMPTY_PATH, false
end

function make_bamparser( bamfile )
   isfile(bamfile) || error("ERROR: --bam parameter used, but cannot find .bam file at $bamfile !")
   isfile(bamfile * ".bai") || error("ERROR: no .bai index found for $bamfile! Cannot set-up random access to $bamfile !")
   println(stderr, "Loading BAM file for random-access: $bamfile")
   bamreadr = open(BAM.Reader, bamfile, index=bamfile * ".bai")
   bamreadr
end

function ispaired_bamfile( bamfile, max_reads=100 )
   parser = open(BAM.Reader, bamfile)
   i = 0
   while !eof(parser) && i < max_reads
      ispaired(first(parser)) && return true
      i += 1
   end
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

function make_collated_parser( filename; flag="-f", force_collated=true )
   check_samtools() || error("Samtools is not properly installed, but is a dependency of BAM input!")
   prefix = splitext(basename(filename))[1]
   head = SAM.Reader(open(`samtools view -H $filename`, "r"))
   iscollated = has_collated_header(head)

   if !iscollated || force_collated
      println(stderr, "Loading BAM file and collating with samtools...")
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
                            mod )
   concordant    = 0
   nonconcordant = 0
   singletons    = 0
   numintrons    = 0

   mean_readlen  = 0.0
  
   leftmate  = eltype(parser)()
   rightmate = eltype(parser)()

   maxdist = 5000000
   # go through all reads, add multi mapping reads
   # after 100K reads calculate ratio and resize multi_reads
   # assuming a 10M library size
   i = 1

   data = AlignData()

   while !eof( parser )

      read!( parser, leftmate )
      read!( parser, rightmate )
      
      while tempname(leftmate) != tempname(rightmate)
         # Process leftmate as singleton
         #align and count
         leftpath, leftvalid = align_bam( lib.coords, data, leftmate )
         if leftvalid
            leftseq = sequence(leftmate) 
            if !isstrandfwd(leftmate)
               BioSequences.reverse_complement!(leftseq)
               reverse!(leftpath)
            end
            biasval = count!( mod, leftseq )
            count!( quant, leftpath, biasval )
            @fastmath mean_readlen += (length(leftseq) - mean_readlen) / i
            i += 1
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

      dist = leftposition(rightmate) - rightposition(leftmate)
      if !isleftmate(leftmate) || 
         !isrightmate(rightmate) ||
         refname(leftmate) != refname(rightmate) ||
         abs(dist) > maxdist ||
         isstrandfwd(leftmate) != isstrandfwd(rightmate)

         # Non-concordant reads
         nonconcordant += 1

      else

         @assert dist > 0
         # Process left/right first_mates
         leftpath, leftvalid   = align_bam( lib.coords, data, leftmate )
         rightpath, rightvalid = align_bam( lib.coords, data, rightmate )

         leftseq  = sequence(leftmate)
         rightseq = sequence(rightmate)

         if isstrandfwd(leftmate)
            biasval = count!( mod, leftseq, rightseq )
         else
            BioSequences.reverse_complement!(leftseq)
            BioSequences.reverse_complement!(rightseq)
            biasval = count!( mod, rightseq, leftseq )
            leftpath, rightpath = reverse(rightpath), reverse(leftpath)
         end
         count!( quant, leftpath, rightpath, biasval )

         @fastmath mean_readlen += (length(leftseq) - mean_readlen) / i
         i += 1
         @fastmath mean_readlen += (length(rightseq) - mean_readlen) / i
         i += 1

         concordant += 1
      end

      if i % 10000000 == 0
         println(i)
      end

      empty!(leftmate)
      empty!(rightmate)
   end
   concordant, nonconcordant, singletons, mean_readlen
end

# for graph.jl:  TODO (not yet implemented for additional RI)
#function is_retained()
# if ! current is intron and usebam=true
    #return false
# end
# get overlapping unspliced read count
# intron_expr = sum / length
# if intron_expr > gene.exon_expr * ratio
    # return true
# else
    # return false
# end
#end



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