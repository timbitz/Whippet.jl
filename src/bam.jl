
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

isskipchar( c::Char ) = c == 'N'
introncount( rec::R ) where R <: Union{SAM.Record, BAM.Record} = count(isskipchar, cigar(rec))
isspliced( rec::R ) where R <: Union{SAM.Record, BAM.Record} = introncount(rec) > 0
isstrandfwd( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x010 == 0
ispaired( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x1 != 0
isleftmate( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x40 != 0
isrightmate( rec::R ) where R <: Union{SAM.Record, BAM.Record} = flag(rec) & 0x80 != 0


function extract_introns( rec::R ) where R <: Union{SAM.Record, BAM.Record}
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


function read_collated_bam( parser; paired_mode=true, niter=100000 )
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

   while !eof( parser ) && i < niter

      read!( parser, leftmate )
      read!( parser, rightmate )
      
      while tempname(leftmate) != tempname(rightmate)
         # Process leftmate as singleton
         #align and count

         singletons += 1
         i += 1

         # Now shift
         leftmate,rightmate = rightmate,leftmate
         empty!(rightmate)
         eof(parser) && break
         read!(parser, rightmate)
         println("shifting singleton at $i")
      end

      absdist = abs(rightposition(leftmate) - leftposition(rightmate))
      if !isleftmate(leftmate) || 
         !isrightmate(rightmate) ||
         refname(leftmate) != refname(rightmate) ||
         absdist > maxdist

         # Non-concordant reads
         println("non-concordant pair at $i")
      end

      # Process left/right first_mates
      if introncount(leftmate) > 0
         introns = extract_introns(leftmate)
         numintrons += length(introns)
      end

      if introncount(rightmate) > 0
         introns = extract_introns(rightmate)
         numintrons += length(introns)
      end

      concordant += 1
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