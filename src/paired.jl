
@inline score{A <: UngappedAlignment}( one::A, two::A  ) = @fastmath score( one ) + score( two )
@inline identity{A <: UngappedAlignment}( one::A, two::A, onelen::Int, twolen::Int ) = @fastmath score( one, two ) / ( onelen + twolen )

@inline function sorted_offsets( sa::UnitRange, index::FMIndexes.FMIndex )
   res = Vector{Int}(length(sa))
   ridx = 1
   for sidx in sa
      res[ridx] = FMIndexes.sa_value( sidx, index ) + 1
      ridx += 1
   end
   sort!(res)
   res
end

# Paired version.
function splice_by_score!{A <: UngappedAlignment}( one::Vector{A}, two::Vector{A}, threshold, buffer )
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

function splice_by_identity!{A <: UngappedAlignment}( one::Vector{A}, two::Vector{A}, threshold, buffer, onelen, twolen )
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

   const fwd_seed,fwd_readloc,strand = seed_locate( p, lib.index, fwd, offset_left=anchor_left, both_strands=!p.is_stranded )
   
   if (p.is_pair_rc && strand) || (!p.is_pair_rc && !strand)
      reverse_complement!(rev)
      const rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=!anchor_left, both_strands=false)
   else
      const rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=anchor_left, both_strands=false)
   end
   const fwd_len = length(fwd.sequence)
   const rev_len = length(rev.sequence)

   if !strand
      reverse_complement!(fwd)
      ispos = false
   end

   fwd_res      = Nullable{Vector{SGAlignment}}()
   rev_res      = Nullable{Vector{SGAlignment}}()
   maxscore = 0.0

   const fwd_sorted = sorted_offsets( fwd_seed, lib.index )
   const rev_sorted = sorted_offsets( rev_seed, lib.index )

   fidx,ridx = 1,1

   while fidx <= length(fwd_sorted) && ridx <= length(rev_sorted)

      const dist = abs( fwd_sorted[fidx] - rev_sorted[ridx] )
      if dist < p.seed_pair_range # inside allowable pair distance

         const geneind = convert( CoordInt, searchsortedlast( lib.offset, fwd_sorted[fidx] ) )
         const next_offset= geneind < length(lib.offset) ? lib.offset[geneind+1] : 
                                                           lib.offset[geneind] + length(lib.graphs[geneind].seq)
         if lib.offset[geneind] <= rev_sorted[ridx] < next_offset

            fwd_aln = _ungapped_align( p, lib, fwd, fwd_sorted[fidx], fwd_readloc; ispos=ispos, geneind=geneind )
            rev_aln = _ungapped_align( p, lib, rev, rev_sorted[ridx], rev_readloc; ispos=!ispos, geneind=geneind )

            @fastmath const scvar = identity( fwd_aln, rev_aln, fwd_len, rev_len )

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

   # if !stranded and no valid alignments, run reverse complement
#=   if ispos && !p.is_stranded && (isnull( fwd_res ) && isnull( rev_res ))
      reverse_complement!( fwd.sequence )
      reverse_complement!( rev.sequence )
      reverse!( fwd.quality )
      reverse!( rev.quality )
      fwd_res,rev_res = ungapped_align( p, lib, fwd, rev, ispos=false, anchor_left=!anchor_left )
   end=#

   fwd_res,rev_res
end


## Extension of count! for paired end counting from quant.jl
function count!( graphq::GraphLibQuant, fwd::SGAlignment, rev::SGAlignment; val=1.0, used=IntSet() )
   (fwd.isvalid == true && rev.isvalid == true) || return
   const init_gene = fwd.path[1].gene
   const rev_gene = rev.path[1].gene

   (init_gene == rev_gene) || return

   sgquant   = graphq.quant[ init_gene ]

   graphq.count[ init_gene ] += 1

   if length(fwd.path) == 1
      # access node -> [ SGNode( gene, *node* ) ]
      sgquant.node[ fwd.path[1].node ] += val
      push!( used, fwd.path[1].node )
   else
      # Otherwise, lets step through pairs of nodes and add val to those edges
      for n in 1:(length(fwd.path)-1)
         # trans-spicing off->
         fwd.path[n].gene != init_gene && continue
         fwd.path[n+1].gene != init_gene && continue
         const lnode = fwd.path[n].node
         const rnode = fwd.path[n+1].node
         push!( used, lnode )
         push!( used, rnode )
         if lnode < rnode
            interv = Interval{ExonInt}( lnode, rnode )
            sgquant.edge[ interv ] = get( sgquant.edge, interv, IntervalValue(0,0,0.0) ).value + val
         elseif lnode >= rnode
            sgquant.circ[ (lnode, rnode) ] = get( sgquant.circ, (lnode,rnode), 0.0) + val
         end
      end
   end

   # now conditionally add rev if not already used.
   if length(rev.path) == 1
      # access node -> [ SGNode( gene, *node* ) ]
      if !(rev.path[1].node in used)
         sgquant.node[ rev.path[1].node ] += val
      end
   else
      # Otherwise, lets step through pairs of nodes and add val to those edges
      for n in 1:(length(rev.path)-1)
         # trans-spicing off->
         rev.path[n].gene != init_gene && continue
         rev.path[n+1].gene != init_gene && continue
         const lnode = rev.path[n].node
         const rnode = rev.path[n+1].node
         (lnode in used && rnode in used) && continue
         if lnode < rnode
            interv = Interval{ExonInt}( lnode, rnode )
            sgquant.edge[ interv ] = get( sgquant.edge, interv, IntervalValue(0,0,0.0) ).value + val
         elseif lnode >= rnode
            sgquant.circ[ (lnode, rnode) ] = get( sgquant.circ, (lnode,rnode), 0.0) + val
         end
      end
   end
   #empty!( used )

end
