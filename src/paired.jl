
function seed_intersect( fwd_sa::UnitRange, rev_sa::UnitRange, index::FMIndexes.FMIndex, max_dist::Int )
   
end

# paired-end seed_locate
function seed_locate( p::AlignParam, index::FMIndex, fwd::SeqRecord, rev::SeqRecord; fwd_offset_left=true, rev_offset_left=false )
   const def_sa = 2:1
   ctry = 1

   if fwd_offset_left
      const fwd_increment = p.seed_inc
      fwd_curpos = p.seed_buffer
   else
      const fwd_increment = -p.seed_inc
      fwd_curpos = length(fwd.seq) - p.seed_length - p.seed_buffer
   end

   if rev_offset_left
      const rev_increment = p.seed_inc
      rev_curpos = p.seed_buffer
   else
      const rev_increment = -p.seed_inc
      rev_curpos = length(rev.seq) - p.seed_length - p.seed_buffer
   end
   
   const fwd_maxright = length(fwd.seq) - p.seed_length - p.seed_buffer
   const rev_maxright = length(rev.seq) - p.seed_length - p.seed_buffer

   while( ctry <= p.seed_try && p.seed_buffer <= fwd_curpos <= fwd_maxright && 
                                p.seed_buffer <= rev_curpos <= rev_maxright )

      const fwd_sa = FMIndexes.sa_range( fwd.seq[fwd_curpos:(fwd_curpos+p.seed_length-1)], index )
      const rev_sa = FMIndexes.sa_range( rev.seq[rev_curpos:(rev_curpos+p.seed_length-1)], index )

      const cnt = length(sa)
      ctry += 1
      #println("$sa, curpos: $curpos, cnt: $cnt, try: $(ctry-1)")
      if cnt == 0 || cnt > p.seed_tolerate
         curpos += increment
      else
         return sa,curpos
      end
   end
   def_sa,curpos
end

function ungapped_align( p::AlignParam, lib::GraphLib, fwd::SeqRecord, rev::SeqRecord; ispos=true, anchor_left=true )

   const fwd_seed,fwd_readloc = seed_locate( p, lib.index, fwd, offset_left=anchor_left )
   
   if p.is_pair_rc
      reverse_complement!(rev.seq)
      reverse!(rev.metadata.quality)
      const rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=false )
   else
      const rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=anchor_left )
   end

   res      = Nullable{Vector{SGAlignment}}()
   maxscore = 0.0
   for sidx in 1:length(seed)
      const s = FMIndexes.sa_value( seed[sidx], lib.index ) + 1
      const geneind = search_sorted( lib.offset, convert(Coordint, s), lower=true )
      align = ungapped_fwd_extend( p, lib, convert(Coordint, geneind),
                                   s - lib.offset[geneind] + p.seed_length,
                                   read, readloc + p.seed_length, ispos=ispos )

      align = ungapped_rev_extend( p, lib, convert(Coordint, geneind),
                                   s - lib.offset[geneind] - 1,
                                   read, readloc - 1, ispos=ispos, align=align,
                                   nodeidx=align.path[1].node )
      if align.isvalid
         if isnull( res )
            res = Nullable(SGAlignment[ align ])
         else
            # new best score
            const scvar = score(align)
            if scvar > maxscore
               # better than threshold for all of the previous scores
               if scvar - maxscore > p.score_range
                  length( res.value ) >= 1 && empty!( res.value )
               else # keep at least one previous score
                  splice_by_score!( res.value, scvar, p.score_range )
               end
               maxscore = scvar
               push!( res.value, align )
            else # new score is lower than previously seen
               if maxscore - scvar <= p.score_range # but tolerable
                  push!( res.value, align )
               end
            end # end score vs maxscore
         end # end isnull 
      end # end isvalid 
   end
   # if !stranded and no valid alignments, run reverse complement
   if ispos && !p.is_stranded && isnull( res )
#      read.seq = Bio.Seq.unsafe_complement!(Bio.Seq.reverse( read.seq ))
      reverse_complement!( read.seq )
      reverse!( read.metadata.quality )
      res = ungapped_align( p, lib, read, ispos=false, anchor_left=!anchor_left )
   end
   res
end
