# requires
# using FMIndexes
# include("graph.jl")
# include("index.jl")

import Bio.Seq.SeqRecord

immutable AlignParam
   mismatches::Int      # Allowable mismatches
   kmer_size::Int       # Minimum number of matches to extend past an edge
   seed_try::Int        # Starting number of seeds
   seed_tolerate::Int   # Allow at most _ hits for a valid seed
   seed_min_qual::Int   # Absolute lowest quality score to allow for worst base in seed
   seed_length::Int     # Seed Size
   seed_buffer::Int     # Ignore first and last _ bases
   seed_inc::Int        # Incrementation for subsequent seed searches
   seed_pair_range::Int # Seeds for paired-end reads must match within _ nt of each other
   score_range::Int     # Scoring range to be considered repetitive
   score_min::Int       # Minimum score for a valid alignment
   is_stranded::Bool    # Is input data + strand only?
   is_paired::Bool      # Paired end data?
   is_pair_rc::Bool     # Is mate 2 the reverse complement of mate 1
   is_trans_ok::Bool    # Do we allow edges from one gene to another
   is_circ_ok::Bool     # Do we allow back edges
end

AlignParam( ispaired = false ) = ispaired ? AlignParam( 2, 9, 2, 4, 4, 18, 5, 18, 2500, 25, 45, 
                                                        false, true, true, false, true ) :
                                            AlignParam( 2, 9, 2, 4, 4, 18, 5, 18, 1000, 10, 45, 
                                                        false, false, true, false, true )

# load from command line args
@inline function AlignParam( args::Dict{AbstractString,Any}, ispaired=false; kmer=9 )
   const aln = AlignParam( args["mismatches"], kmer,
                           args["seed-try"],
                           args["seed-tol"], 4,
                           args["seed-len"],
                           args["seed-buf"],
                           args["seed-inc"],
                           args["pair-range"], 10,
                           args["score-min"],
                           args["stranded"], ispaired,
                           args["rev-pair"], false,
                           args["no-circ"] )
   aln
end

abstract UngappedAlignment

type SGAlignment <: UngappedAlignment
   matches::Int
   mismatches::Float64
   offset::Int
   path::Vector{SGNode}
   strand::Bool
   isvalid::Bool
end

SGAlignment() = SGAlignment(0, 0, 0, SGNode[], true, false)

typealias SGAlignVec Nullable{Vector{SGAlignment}}

const DEF_ALIGN = SGAlignment(0, 0, 0, SGNode[], true, false)

@inline score{A <: UngappedAlignment}( align::A ) = @fastmath align.matches - align.mismatches 

Base.(:>)( a::SGAlignment, b::SGAlignment ) = >( score(a), score(b) )
Base.(:<)( a::SGAlignment, b::SGAlignment ) = <( score(a), score(b) )
Base.(:(>=))( a::SGAlignment, b::SGAlignment ) = >=( score(a), score(b) )
Base.(:(<=))( a::SGAlignment, b::SGAlignment ) = <=( score(a), score(b) )

# add prob of being accurate base to mismatch, rather than integer.
@inline phred_to_prob( phred::Int8 ) = @fastmath 1-10^(-phred/10)

function Base.empty!( align::SGAlignment ) 
   align.matches = 0 
   align.mistmatches = 0 
   align.offset = 0 
   empty!( align.path ) 
   align.strand = true 
   align.isvalid = false 
end

@inline function seed_locate( p::AlignParam, index::FMIndex, read::SeqRecord; offset_left=true, both_strands=true )
   const def_sa = 2:1
   ctry = 1
   ispos = true
   if offset_left
      const increment = p.seed_inc
      const increment_sm = 1
      curpos = p.seed_buffer
   else
      const increment = -p.seed_inc
      const increment_sm = -1
      curpos = length(read.seq) - p.seed_length - p.seed_buffer
   end
   const maxright = length(read.seq) - p.seed_length - p.seed_buffer
   while( ctry <= p.seed_try && p.seed_buffer <= curpos <= maxright )
      if minimum( read.metadata.quality[curpos:(curpos+p.seed_length-1)] ) <= p.seed_min_qual
         curpos += increment_sm
         continue
      end
      const curseq = read.seq[curpos:(curpos+p.seed_length-1)]
      const sa = FMIndexes.sa_range( curseq, index )
      const cnt = length(sa)
      ctry += 1
      #println("$sa, curpos: $curpos, cnt: $cnt, try: $(ctry-1)")
      if cnt == 0 || cnt > p.seed_tolerate
         if both_strands
            reverse_complement!(curseq)
            const rev_sa = FMIndexes.sa_range( curseq, index )
            const rev_cnt = length(rev_sa)
            if 0 < rev_cnt <= p.seed_tolerate
               const rev_pos = length(read.seq)-(curpos+p.seed_length-1)+1
               return rev_sa,rev_pos,false
            end
         end
      else
         return sa,curpos,true
      end
      curpos += increment
   end
   def_sa,curpos,true
end

function splice_by_score!{A <: UngappedAlignment}( arr::Vector{A}, threshold, buffer )
   i = 1
   while i <= length( arr )
      if threshold - score( arr[i] ) > buffer
         splice!( arr, i )
         i -= 1
      end
      i += 1
   end
end

@inline function _ungapped_align( p::AlignParam, lib::GraphLib, read::SeqRecord, indx::Int, readloc::Int;
                                  ispos=true, geneind=convert(Coordint, searchsortedlast( lib.offset, indx )) )

   align = ungapped_fwd_extend( p, lib, geneind,
                                indx - lib.offset[geneind] + p.seed_length,
                                read, readloc + p.seed_length, ispos=ispos )

   align = ungapped_rev_extend( p, lib, geneind,
                                indx - lib.offset[geneind] - 1,
                                read, readloc - 1, ispos=ispos, align=align,
                                nodeidx=align.path[1].node )

   align
end

function ungapped_align( p::AlignParam, lib::GraphLib, read::SeqRecord; ispos=true, anchor_left=true )
   const seed,readloc,pos = seed_locate( p, lib.index, read, offset_left=anchor_left, both_strands=!p.is_stranded )
   #@bp
   if !pos
      reverse_complement!( read.seq )
      reverse!( read.metadata.quality )
      ispos = false
   end

   res      = Nullable{Vector{SGAlignment}}()
   maxscore = 0.0
#   for s in FMIndexes.LocationIterator( seed, lib.index ) # rm extra malloc
   for sidx in 1:length(seed)
      const s = FMIndexes.sa_value( seed[sidx], lib.index ) + 1
      #println("$(read.seq[readloc:(readloc+75)])\n$(lib.graphs[geneind].seq[(s-lib.offset[geneind]):(s-lib.offset[geneind])+50])")
      #@bp

      const align = _ungapped_align( p, lib, read, s, readloc, ispos=ispos )

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
   #=if ispos && !p.is_stranded && isnull( res )
#      read.seq = Bio.Seq.unsafe_complement!(Bio.Seq.reverse( read.seq ))
      reverse_complement!( read.seq )
      reverse!( read.metadata.quality )
      res = ungapped_align( p, lib, read, ispos=false, anchor_left=!anchor_left )
   end =#
   res
end


# This is the main ungapped alignment extension function in the --> direction
# Returns: SGAlignment
function ungapped_fwd_extend( p::AlignParam, lib::GraphLib, geneind::Coordint, sgidx::Int, 
                                read::SeqRecord, ridx::Int; ispos=true,
                                align::SGAlignment=SGAlignment(p.seed_length,0,sgidx,SGNode[],ispos,false),
                                nodeidx=searchsortedlast(lib.graphs[geneind].nodeoffset,Coordint(sgidx)) )
   const sg       = lib.graphs[geneind]
   const readlen  = length(read.seq)

   #curedge  = nodeidx

   passed_edges = Nullable{Vector{UInt32}}() # don't allocate array unless needed
   passed_extend = 0
   passed_mismat = 0.0

   push!( align.path, SGNode( geneind, nodeidx ) ) # starting node

   while( align.mismatches <= p.mismatches && ridx <= readlen && sgidx <= length(sg.seq) )
      if read.seq[ridx] == sg.seq[sgidx]
         # match
         align.matches += 1
         passed_extend  += 1
      elseif (UInt8(sg.seq[sgidx]) & 0b100) == 0b100 && !(sg.seq[sgidx] == SG_N) # L,R,S
         const curedge::UInt32 = nodeidx+1
         #@bp
         if     sg.edgetype[curedge] == EDGETYPE_LR &&
                sg.nodelen[curedge]  >= p.kmer_size && # 'LR' && nodelen >= K
                readlen - ridx + 1   >= p.kmer_size
               # check edgeright[curnode+1] == read[ nextkmer ]
               # move forward K, continue 
               # This is a std exon-exon junction, if the read is an inclusion read
               # then we simply jump ahead, otherwise lets try to find a compatible right
               # node
               const rseq = read.seq[ridx:(ridx+p.kmer_size-1)]
               Bio.Seq.hasn( rseq ) && break
               const rkmer = DNAKmer{p.kmer_size}( rseq )
               if sg.edgeright[curedge] == rkmer
                  align.matches += p.kmer_size
                  ridx    += p.kmer_size - 1
                  sgidx   += 1 + p.kmer_size
                  nodeidx += 1
                  push!( align.path, SGNode( geneind, nodeidx ) )
                  align.isvalid = true
               else
                  align = spliced_fwd_extend( p, lib, geneind, curedge, read, ridx, rkmer, align )
                  break
               end
         elseif sg.edgetype[curedge] == EDGETYPE_LR || 
                sg.edgetype[curedge] == EDGETYPE_LL # 'LR' || 'LL'
               # if length(edges) > 0, push! edgematch then set to 0
               # push! edge,  (here we try to extend fwd, but keep track
               # of potential edges we pass along the way.  When the alignment
               # dies if we have not had sufficient matches we then explore
               # those edges. 
               if isnull(passed_edges) # now we have to malloc
                  passed_edges = Nullable(Vector{UInt32}())
               end
               passed_extend = 0
               passed_mismat = 0
               push!(get(passed_edges), curedge)
               #@bp
               sgidx   += 1
               ridx    -= 1
               nodeidx += 1
               push!( align.path, SGNode( geneind, nodeidx ) )
         elseif sg.edgetype[curedge] == EDGETYPE_LS && # 'LS'
                readlen - ridx + 1   >= p.kmer_size
               # obligate spliced_extension
               const rseq = read.seq[ridx:(ridx+p.kmer_size-1)]
               Bio.Seq.hasn( rseq ) && break
               const rkmer = DNAKmer{p.kmer_size}( rseq )
               align = spliced_fwd_extend( p, lib, geneind, curedge, read, ridx, rkmer, align )
               break
         elseif sg.edgetype[curedge] == EDGETYPE_SR #||
                #@bp
               # end of alignment
               break # ?
         else #'RR' || 'SL' || 'RS'
               # ignore 'RR' and 'SL'
               sgidx += 1
               ridx  -= 1 # offset the lower ridx += 1
               nodeidx += 1
               nodeidx <= length(sg.nodelen) && push!( align.path, SGNode( geneind, nodeidx ) )  
               #@bp
         end
         # ignore 'N'
      else 
         # mismatch
         const prob  = phred_to_prob( read.metadata.quality[ridx] )
         @fastmath align.mismatches += prob
         @fastmath passed_mismat += prob
      end
      ridx  += 1
      sgidx += 1
      #print(" $(read.seq[ridx-1]),$ridx\_$(sg.seq[sgidx-1]),$sgidx ")
   end


   # if edgemat < K, spliced_extension for each in length(edges)
   if !isnull(passed_edges)
      if passed_extend < p.kmer_size
         # go back.
         const ext_len = Int[sg.nodelen[ i ] for i in get(passed_edges)]
         ext_len[end] = passed_extend
         rev_cumarray!(ext_len)
         #@bp
         for c in length(get(passed_edges)):-1:1  #most recent edge first
            if ext_len[c] >= p.kmer_size
               align.isvalid = true
               break
            else
               pop!( align.path ) # extension into current node failed
               @fastmath align.mismatches -= passed_mismat
               align.matches -= ext_len[c]
               const cur_ridx = ridx - (sgidx - sg.nodeoffset[ get(passed_edges)[c] ])
               (cur_ridx + p.kmer_size - 1) <= readlen || continue
            
               #lkmer = DNAKmer{p.kmer_size}(read.seq[(ridx-p.kmer_size):(ridx-1)])
               const rseq  = read.seq[cur_ridx:(cur_ridx+p.kmer_size-1)]
               Bio.Seq.hasn( rseq ) && continue
               const rkmer = DNAKmer{p.kmer_size}( rseq ) 
               res_align = spliced_fwd_extend( p, lib, geneind, get(passed_edges)[c], 
                                               read, cur_ridx, rkmer, align )
                  #println("$res_align")
               if length(res_align.path) > length(align.path) # we added a valid node
                  align = res_align > align ? res_align : align
               end
            end
         end
      else
         align.isvalid = true
      end
   end

   #TODO alignments that hit the end of the read but < p.kmer_size beyond an
   # exon-exon junction should still be valid no?  Even if they aren't counted
   # as junction spanning reads.  Currently they would not be.
   if score(align) >= p.score_min && length(align.path) == 1  
      align.isvalid = true 
   end
   align
end



function spliced_fwd_extend{T,K}( p::AlignParam, lib::GraphLib, geneind::Coordint, edgeind::UInt32, 
                                  read::SeqRecord, ridx::Int, rkmer::Bio.Seq.Kmer{T,K}, align::SGAlignment )
   # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
   # returns right nodes with matching genes
   isdefined( lib.edges.right, kmer_index(rkmer) ) || return align

   best = align
   const left_kmer_ind = kmer_index(lib.graphs[geneind].edgeleft[edgeind])
   isdefined( lib.edges.left, left_kmer_ind ) || return align
   const right_nodes = lib.edges.left[ left_kmer_ind ] ∩ lib.edges.right[ kmer_index(rkmer) ]

   for rn in right_nodes
      rn.gene == align.path[1].gene || continue
      const rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node]
      res_align::SGAlignment = ungapped_fwd_extend( p, lib, rn.gene, Int(rn_offset), read, ridx, 
                                                    align=deepcopy(align), nodeidx=rn.node )
      best = res_align > best ? res_align : best
   end
   if score(best) >= score(align) + p.kmer_size
      best.isvalid = true
   end
   best
end

# Function Cumulative Array
@inline function rev_cumarray!{T<:Vector}( array::T )
   for i in (length(array)-1):-1:1
      array[i] += array[i+1]
   end 
end

# This is the main ungapped alignment extension function in the <-- direction
# Returns: SGAlignment
function ungapped_rev_extend( p::AlignParam, lib::GraphLib, geneind::Coordint, sgidx::Int, 
                                read::SeqRecord, ridx::Int; ispos=true,
                                align::SGAlignment=SGAlignment(p.seed_length,0,sgidx,SGNode[],ispos,false),
                                nodeidx=searchsortedlast(lib.graphs[geneind].nodeoffset,Coordint(sgidx)) )
   const sg       = lib.graphs[geneind]
   const readlen  = length(read.seq)

   #leftnode  = nodeidx

   passed_edges = Nullable{Vector{UInt32}}() # don't allocate array unless needed
   passed_extend = 0
   passed_mismat = 0.0

   if align.path[1] != SGNode( geneind, nodeidx )
      unshift!( align.path, SGNode( geneind, nodeidx ) ) # starting node if not already there
   end

   while( align.mismatches <= p.mismatches && ridx > 0 && sgidx > 0 )
      #@bp
      if read.seq[ridx] == sg.seq[sgidx]
         # match
         align.matches += 1
         passed_extend  += 1
      elseif (UInt8(sg.seq[sgidx]) & 0b100) == 0b100 && !(sg.seq[sgidx] == SG_N) # L,R,S
         const leftnode = nodeidx - 1
         #@bp
         if     sg.edgetype[nodeidx] == EDGETYPE_LR &&
                sg.nodelen[leftnode] >= p.kmer_size && # 'LR' && nodelen >= K
                                ridx >= p.kmer_size 
                
               const rseq  = read.seq[(ridx-p.kmer_size+1):ridx]
               Bio.Seq.hasn( rseq ) && break 
               const lkmer = DNAKmer{p.kmer_size}( rseq )
               if sg.edgeleft[nodeidx] == lkmer
                  align.matches += p.kmer_size
                  ridx    -= p.kmer_size - 1 
                  sgidx   -= p.kmer_size + 1 
                  nodeidx -= 1
                  unshift!( align.path, SGNode( geneind, nodeidx ) )
                  align.isvalid = true
               else
                  align = spliced_rev_extend( p, lib, geneind, convert(UInt32,nodeidx), read, ridx, lkmer, align ) 
                  break
               end
         elseif sg.edgetype[nodeidx] == EDGETYPE_LR || 
                sg.edgetype[nodeidx] == EDGETYPE_RR # 'LR' || 'RR'
                
               if isnull(passed_edges) # now we have to malloc
                  passed_edges = Nullable(Vector{UInt32}())
               end
               passed_extend = 0
               passed_mismat = 0
               push!(get(passed_edges), nodeidx)
              # @bp
               sgidx   -= 1
               ridx    += 1
               nodeidx -= 1
               unshift!( align.path, SGNode( geneind, nodeidx ) )
         elseif sg.edgetype[nodeidx] == EDGETYPE_SR && # 'SR'
                                ridx >= p.kmer_size
               # obligate spliced_extension
               const rseq  = read.seq[(ridx-p.kmer_size+1):ridx]
               Bio.Seq.hasn( rseq ) && break
               const lkmer = DNAKmer{p.kmer_size}( rseq )
               align = spliced_rev_extend( p, lib, geneind, convert(UInt32,nodeidx), read, ridx, lkmer, align )
               break
         elseif sg.edgetype[nodeidx] == EDGETYPE_LS #||
               # end of alignment
               break # ?
         else #'LL' || 'RS' || 'SL'
               # ignore 'LL' and 'RS'
               sgidx -= 1
               ridx  += 1 # offset the lower ridx += 1
               nodeidx -= 1
               nodeidx > 0 && unshift!( align.path, SGNode( geneind, nodeidx ) )  
         end
         # ignore 'N'
      else 
         # mismatch
         const prob  = phred_to_prob( read.metadata.quality[ridx] )
         @fastmath align.mismatches += prob
         @fastmath passed_mismat += prob
      end
      ridx  -= 1
      sgidx -= 1
      #print(" $(read.seq[ridx+1]),$ridx\_$(sg.seq[sgidx+1]),$sgidx ")
   end

   # if passed_extend < K, spliced_extension for each in length(edges)
   if !isnull(passed_edges)
      if passed_extend < p.kmer_size
         # go back.
         const ext_len = Int[sg.nodelen[ i ] for i in get(passed_edges)]
         ext_len[end] = passed_extend
         rev_cumarray!(ext_len)
         #@bp
         for c in length(get(passed_edges)):-1:1  #most recent edge first
            if ext_len[c] >= p.kmer_size
               align.isvalid = true
               break
            else
               shift!( align.path ) # extension into current node failed
               @fastmath align.mismatches -= passed_mismat
               align.matches -= ext_len[c]
               const cur_ridx = ridx + (sg.nodeoffset[ get(passed_edges)[c] ] - sgidx - 2)
               (cur_ridx - p.kmer_size) > 0 || continue
            
               const rseq  = read.seq[(cur_ridx-p.kmer_size):(cur_ridx-1)]
               Bio.Seq.hasn( rseq ) && continue
               const lkmer = DNAKmer{p.kmer_size}( rseq )
               res_align = spliced_rev_extend( p, lib, geneind, get(passed_edges)[c], 
                                                  read, cur_ridx-1, lkmer, align )
               if length(res_align.path) > length(align.path) # we added a valid node
                  align = res_align > align ? res_align : align
               end
            end
         end
      else
         align.isvalid = true
      end
   end

   #TODO alignments that hit the end of the read but < p.kmer_size beyond an
   # exon-exon junction should still be valid no?  Even if they aren't counted
   # as junction spanning reads.  Currently they would not be.
   if score(align) >= p.score_min && length(align.path) == 1  
      align.isvalid = true 
   end
   if sgidx < align.offset
      align.offset = sgidx
   end
   align
end

function spliced_rev_extend{T,K}( p::AlignParam, lib::GraphLib, geneind::Coordint, edgeind::Coordint,
                                  read::SeqRecord, ridx::Int, lkmer::Bio.Seq.Kmer{T,K}, align::SGAlignment )
   # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
   # returns right nodes with matching genes
   isdefined( lib.edges.left, kmer_index(lkmer) ) || return align

   best = align
   const right_kmer_ind = kmer_index(lib.graphs[geneind].edgeright[edgeind])
   isdefined( lib.edges.right, right_kmer_ind ) || return align
   const left_nodes = lib.edges.right[ right_kmer_ind ] ∩ lib.edges.left[ kmer_index(lkmer) ]

   #println("$left_nodes, $lkmer")
   #@bp

   # do a test for trans-splicing, and reset left_nodes
   for rn in left_nodes
      rn.gene == align.path[1].gene || continue
      const rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node-1] + lib.graphs[rn.gene].nodelen[rn.node-1] - 1
      #println("$(lib.graphs[rn.gene].seq[rn_offset-50:rn_offset])") 
      res_align::SGAlignment = ungapped_rev_extend( p, lib, rn.gene, Int(rn_offset), read, ridx, 
                                                    align=deepcopy(align), nodeidx=rn.node-1 )
      best = res_align > best ? res_align : best
   end
   if score(best) >= score(align) + p.kmer_size
      best.isvalid = true
   end
   best
end
