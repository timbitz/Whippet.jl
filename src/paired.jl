
@inline score( one::A, two::A  ) where {A <: UngappedAlignment} = @fastmath score( one ) + score( two )
@inline identity( one::A, two::A, onelen::Int, twolen::Int ) where {A <: UngappedAlignment} = @fastmath score( one, two ) / ( onelen + twolen )

@inline function sorted_offsets( sa::UnitRange, index::FMIndexes.FMIndex )
   res = Vector{Int}(undef, length(sa))
   ridx = 1
   for sidx in sa
      res[ridx] = FMIndexes.sa_value( sidx, index ) + 1
      ridx += 1
   end
   sort!(res)
   res
end

# Paired version.
function splice_by_score!( one::Vector{A}, two::Vector{A}, threshold, buffer ) where A <: UngappedAlignment
   i = 1
   while i <= length( one )
      if threshold - score( one[i], two[i] ) > buffer
         splice!( one, i )
         splice!( two, i )
         i -= 1
      end
      i += 1
   end
end

function splice_by_identity!( one::Vector{A}, two::Vector{A}, threshold, buffer, onelen, twolen ) where A <: UngappedAlignment
   i = 1
   while i <= length( one )
      if threshold - identity( one[i], two[i], onelen, twolen ) > buffer
         splice!( one, i )
         splice!( two, i )
         i -= 1
      end
      i += 1
   end
end

# Paired ungapped_align
function ungapped_align( p::AlignParam, lib::GraphLib, fwd::FASTQRecord, rev::FASTQRecord; ispos=true, anchor_left=true )

   fwd_seed,fwd_readloc,strand = seed_locate( p, lib.index, fwd, offset_left=anchor_left, both_strands=!p.is_stranded )

   if (p.is_pair_rc && strand) || (!p.is_pair_rc && !strand)
      reverse_complement!(rev)
      rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=!anchor_left, both_strands=false)
   else
      rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=anchor_left, both_strands=false)
   end
   fwd_len = length(fwd.sequence)
   rev_len = length(rev.sequence)

   if !strand
      reverse_complement!(fwd)
      ispos = false
   end

   fwd_res      = Nullable{Vector{SGAlignment}}()
   rev_res      = Nullable{Vector{SGAlignment}}()
   maxscore = 0.0

   fwd_sorted = sorted_offsets( fwd_seed, lib.index )
   rev_sorted = sorted_offsets( rev_seed, lib.index )

   fidx,ridx = 1,1

   while fidx <= length(fwd_sorted) && ridx <= length(rev_sorted)

      dist = abs( fwd_sorted[fidx] - rev_sorted[ridx] )
      if dist < p.seed_pair_range # inside allowable pair distance

         geneind = convert( CoordInt, searchsortedlast( lib.offset, fwd_sorted[fidx] ) )
         next_offset = geneind < length(lib.offset) ? lib.offset[geneind+1] :
                                                           lib.offset[geneind] + length(lib.graphs[geneind].seq)
         if lib.offset[geneind] <= rev_sorted[ridx] < next_offset

            fwd_aln = _ungapped_align( p, lib, fwd, fwd_sorted[fidx], fwd_readloc; ispos=ispos, geneind=geneind )
            rev_aln = _ungapped_align( p, lib, rev, rev_sorted[ridx], rev_readloc; ispos=!ispos, geneind=geneind )

            @fastmath scvar = identity( fwd_aln, rev_aln, fwd_len, rev_len )

            if isvalid(fwd_aln) && isvalid(rev_aln)
               if isnull( fwd_res ) || isnull( rev_res )
                  fwd_res = Nullable(SGAlignment[ fwd_aln ])
                  rev_res = Nullable(SGAlignment[ rev_aln ])
                  maxscore = scvar
               else
               # new best score
                  if scvar > maxscore
                     # better than threshold for all of the previous scores
                     if scvar - maxscore > p.score_range
                        if length( fwd_res.value ) >= 1
                           empty!( fwd_res.value )
                           empty!( rev_res.value )
                        end
                     else # keep at least one previous score
                        splice_by_identity!( fwd_res.value, rev_res.value, scvar, p.score_range, fwd_len, rev_len )
                     end
                     maxscore = scvar
                     push!( fwd_res.value, fwd_aln )
                     push!( rev_res.value, rev_aln )
                  else # new score is lower than previously seen
                     if maxscore - scvar <= p.score_range # but tolerable
                        push!( fwd_res.value, fwd_aln )
                        push!( rev_res.value, rev_aln )
                     end
                  end # end score vs maxscore
               end # end isnull
            end # end isvalid
            ridx += 1
            fidx += 1
            continue
         end # is reverse in same gene as fwd
      end # end dist vs. seed_pair_range
      if fwd_sorted[fidx] > rev_sorted[ridx]
         ridx += 1
      else
         fidx += 1
      end
   end #end while

   fwd_res,rev_res
end
