

@inline function sorted_offsets( sa::UnitRange, index::FMIndexes.FMIndex )
   res = Vector{Int}(length(sa))
   ridx = 1
   for sidx in sa
      res[ridx] = FMIndexes.sa_value( seed[sidx], lib.index ) + 1
      ridx += 1
   end
   sort!(res)
   res
end

function ungapped_align( p::AlignParam, lib::GraphLib, fwd::SeqRecord, rev::SeqRecord; ispos=true, anchor_left=true )

   const fwd_seed,fwd_readloc = seed_locate( p, lib.index, fwd, offset_left=anchor_left )
   
   if p.is_pair_rc
      reverse_complement!(rev.seq)
      reverse!(rev.metadata.quality)
      const rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=!anchor_left )
   else
      const rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=anchor_left )
   end

   fwd_res      = Nullable{Vector{SGAlignment}}()
   rev_res      = Nullable{Vector{SGAlignment}}()
   maxscore = 0.0

   const fwd_sorted = sorted_offsets( fwd_seed, lib.index )
   const rev_sorted = sorted_offsets( rev_seed, lib.index )

   fidx,ridx = 1,1

   while fidx <= length(fwd_sorted) && ridx <= length(rev_sorted)

      const dist = abs( fwd_sorted[fidx] - rev_sorted[ridx] )
      if dist > p.seed_pair_range # outside allowable pair distance
         if fwd_sorted[fidx] > rev_sorted[ridx]
            ridx += 1
         else
            fidx += 1
         end
         continue
      else
         fwd_res = _ungapped_align( p, lib, fwd_res, fwd, fwd_sorted[fidx] )
         rev_res = _ungapped_align( p, lib, rev_res, rev, rev_sorted[ridx] )
      end

   end
   # if !stranded and no valid alignments, run reverse complement
   if ispos && !p.is_stranded && (isnull( fwd_seed ) || isnull( rev_seed ))
      reverse_complement!( fwd.seq )
      reverse_complement!( rev.seq )
      reverse!( fwd.metadata.quality )
      reverse!( rev.metadata.quality )
      fwd_res,rev_res = ungapped_align( p, lib, fwd, rev, ispos=false, anchor_left=!anchor_left )
   end

   fwd_res,rev_res
end

@inline function _ungapped_align( p::AlignParam, lib::GraphLib, read::SeqRecord, indx::Int,
                                     main_res::Nullable{Vector{SGAlignment}}, 
                                  partner_res::Nullable{Vector{SGAlignment}}, maxscore::Int; test_vals=true )

   const geneind = search_sorted( lib.offset, convert(Coordint, indx), lower=true )
   align = ungapped_fwd_extend( p, lib, convert(Coordint, geneind),
                                indx - lib.offset[geneind] + p.seed_length,
                                read, readloc + p.seed_length, ispos=ispos )

   align = ungapped_rev_extend( p, lib, convert(Coordint, geneind),
                                indx - lib.offset[geneind] - 1,
                                read, readloc - 1, ispos=ispos, align=align,
                                nodeidx=align.path[1].node )

   if align.isvalid
      if isnull( main_res )
         main_res = Nullable(SGAlignment[ align ])
      else
            # new best score
         const scvar = score(align)
         if scvar > maxscore
            # better than threshold for all of the previous scores
            if scvar - maxscore > p.score_range
               length( main_res.value ) >= 1 && empty!( main_res.value )
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
   main_res,partner_res
end
