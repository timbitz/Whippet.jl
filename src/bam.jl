
isspliced( rec::BAM.Record ) = occursin(r"N", BAM.cigar(rec))
strandpos( rec::BAM.Record ) = BAM.flag(rec) & 0x010 == 0

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
         if isspliced(rec) && strand == strandpos(rec)
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
   refpos = BAM.position(rec)
   for i in 1:length(op)
      if ismatchop(op[i])
         refpos += len[i]
      elseif isdeleteop(op[i])
         donor = CoordInt(refpos - 1)
         refpos += len[i]
         accep = CoordInt(refpos)
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
