
isspliced( rec::BAM.Record ) = ismatch(r"N", BAM.cigar(i))
strandpos( rec::BAM.Record ) = BAM.flag(rec) & 0x010 == 0

function process_records!( reader::BAM.Reader, seqname::String, range::UnitRange{Int64}, strand::Bool, 
                           exons::CoordTree, novelacc::Dict{K,V}, noveldon::Dict{K,V} ) where {K,V}
   for rec in eachoverlap( reader, seqname, range )
      if isspliced(rec) && strand == strandpos(rec)
         process_spliced_record!( novelacc, noveldon, rec )
      end
   end
end

function process_spliced_record!( novelacc, noveldon, rec::BAM.Record )
   op,len = BAM.cigar_rle( rec )
   refpos = position(rec) - 1
   for i in 1:length(op)
      if ismatchop(op[i])
         refpos += len[i]
      elseif isdeleteop(op[i])
         donor = refpos       #TODO check these
         accep = donor + len[i]
         increment!( noveldon, donor, 1 )
         increment!( novelacc, accep, 1 )
      end
   end
end

# Add to genes type:
# noveldon::Dict{CoordInt,Int}
# novelacc::Dict{CoordInt,Int}
# exon_expr::Float64

#for all introns (refset.jl)
# intron procedure.
# add new splice site(s) to dict
#end

#for all new splice sites (refset.jl)
# delete if splice site is below threshold
# add it to the don/acc lists
#end

#for all exons (in refset.jl)
# add overlapping read count
# add length
#end
#exon_exp = over_read_count / length

# for graph.jl:
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

