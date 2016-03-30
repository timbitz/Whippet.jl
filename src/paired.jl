
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

      const sa = FMIndexes.sa_range( fwd.seq[fwd_curpos:(fwd_curpos+p.seed_length-1)], index )
      

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


