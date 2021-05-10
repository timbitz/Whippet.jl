
const AlignBlocks = Vector{Tuple{CoordInt,CoordInt}}
const AlignNodes  = Vector{GenomicFeatures.Interval{SGNodeIsExon}}

const emptypath   = Vector{SGNodeIsExon}()

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
                             j_a::J, j_b::J) where {I <: Integer, J <: Integer}
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
   length(data.overlaps) > 0 || return emptypath, 0

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
   if all_introns_valid( data, gene )
      return path, true
   end
   empty_path, false
end

function make_bamparser( bamfile )
   isfile(bamfile) || error("ERROR: --bam parameter used, but cannot find .bam file at $bamfile !")
   isfile(bamfile * ".bai") || error("ERROR: no .bai index found for $bamfile! Cannot set-up random access to $bamfile !")
   println(stderr, "Loading BAM file for random-access: $bamfile")
   bamreadr = open(BAM.Reader, bamfile, index=bamfile * ".bai")
   bamreadr
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

function make_collated_parser( filename; flag="-f", force_collated=true )
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
   parser, iscollated
end


function read_collated_bam( parser, lib; paired_mode=true, niter=100000 )
   concordant = 0
   singletons = 0
   numintrons = 0

   leftmate  = eltype(parser)()
   rightmate = eltype(parser)()

   maxdist = 1000000
   # go through all reads, add multi mapping reads
   # after 100K reads calculate ratio and resize multi_reads
   # assuming a 10M library size
   i = 0

   data = AlignData()

   while !eof( parser ) && i < niter

      read!( parser, leftmate )
      read!( parser, rightmate )
      
      while tempname(leftmate) != tempname(rightmate)
         # Process leftmate as singleton
         #align and count
         leftpath, leftvalid = align_bam( lib.coords, data, leftmate )

         singletons += 1
         i += 1

         # Now shift
         leftmate,rightmate = rightmate,leftmate
         empty!(rightmate)
         eof(parser) && break
         read!(parser, rightmate)
         #println("shifting singleton at $i")
      end

      absdist = abs(rightposition(leftmate) - leftposition(rightmate))
      if !isleftmate(leftmate) || 
         !isrightmate(rightmate) ||
         refname(leftmate) != refname(rightmate) ||
         absdist > maxdist

         # Non-concordant reads
         #println("non-concordant pair at $i")
      else

      # Process left/right first_mates
         leftpath, leftvalid   = align_bam( lib.coords, data, leftmate )
         rightpath, rightvalid = align_bam( lib.coords, data, rightmate )
      
         #=if introncount(leftmate) > 0
            introns = extract_introns(leftmate)
            numintrons += length(introns)
         end

         if introncount(rightmate) > 0
            introns = extract_introns(rightmate)
            numintrons += length(introns)
         end=#

         concordant += 1
      end

      i += 1
      if i % 10000000 == 0
         println(i)
      end

      empty!(leftmate)
      empty!(rightmate)
   end
   concordant, singletons, numintrons
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



function process_records!( reader::BAM.Reader, seqname::String, range::UnitRange{Int64}, strand::Bool,
                           exons::CoordTree, known, oneknown::Bool,
                           novelacc::Dict{K,V}, noveldon::Dict{K,V} ) where {K,V}
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

function process_spliced_record!( novelacc::Dict{K,V}, noveldon::Dict{K,V}, rec::BAM.Record,
                                  known, oneknown::Bool ) where {K,V}
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