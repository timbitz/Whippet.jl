
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
   score_range::Float64 # Identity scoring range to be considered repetitive
   score_min::Float64   # Minimum percent identity score for a valid alignment
   is_stranded::Bool    # Is input data + strand only?
   is_paired::Bool      # Paired end data?
   is_pair_rc::Bool     # Is mate 2 the reverse complement of mate 1
   is_trans_ok::Bool    # Do we allow edges from one gene to another
   is_circ_ok::Bool     # Do we allow back edges
end

AlignParam( ispaired = false ) = ispaired ? AlignParam( 3, 9, 2, 4, 4, 18, 5, 18, 2500, 0.05, 0.7, 
                                                        false, true, true, false, true ) :
                                            AlignParam( 3, 9, 2, 4, 4, 18, 5, 18, 1000, 0.05, 0.7, 
                                                        false, false, true, false, true )

# load from command line args
@inline function AlignParam( args::Dict{AbstractString,Any}, ispaired=false; kmer=9 )
   const aln = AlignParam( args["mismatches"], kmer,
                           args["seed-try"],
                           args["seed-tol"], 4,
                           args["seed-len"],
                           args["seed-buf"],
                           args["seed-inc"],
                           args["pair-range"], 0.05,
                           args["score-min"],
                           args["stranded"], ispaired,
                          !args["pair-same-strand"], false,
                           args["circ"] )
   aln
end

abstract UngappedAlignment

type SGAlignment <: UngappedAlignment
   matches::UInt16
   mismatches::Float32
   offset::UInt32
   offsetread::UInt16
   path::Vector{SGNode}
   strand::Bool
   isvalid::Bool
end

SGAlignment() = SGAlignment(0x0000, zero(Float32), zero(UInt32), zero(UInt16), SGNode[], true, false)

typealias SGAlignVec Nullable{Vector{SGAlignment}}

const DEF_ALIGN = SGAlignment(0x0000, zero(Float32), zero(UInt32), zero(UInt16), SGNode[], true, false)

@inline score{A <: UngappedAlignment}( align::A ) = @fastmath align.matches - align.mismatches
@inline identity{A <: UngappedAlignment, I <: Integer}( align::A, readlen::I ) = convert(Float32, @fastmath ( align.matches - align.mismatches ) / readlen)

Base.:>( a::SGAlignment, b::SGAlignment ) = >( score(a), score(b) )
Base.:<( a::SGAlignment, b::SGAlignment ) = <( score(a), score(b) )
Base.:>=( a::SGAlignment, b::SGAlignment ) = >=( score(a), score(b) )
Base.:<=( a::SGAlignment, b::SGAlignment ) = <=( score(a), score(b) )
Base.isless( a::SGAlignment, b::SGAlignment ) = a < b

# add prob of being accurate base to mismatch, rather than integer.
@inline phred_to_prob( phred::Int8 ) = convert(Float32, @fastmath 1-10^(-phred/10))

function Base.empty!( align::SGAlignment ) 
   align.matches = 0x0000
   align.mismatches = zero(Float32)
   align.offset = zero(UInt32)
   align.offsetread = zero(UInt16)
   empty!( align.path ) 
   align.strand = true 
   align.isvalid = false 
end

@inline function seed_locate( p::AlignParam, index::FMIndex, read::SeqRecord; offset_left::Bool=true, both_strands::Bool=true )
   const def_sa = 2:1
   const readlen = convert(UInt16, length(read.seq))
   ctry = 1
   #ispos = true
   if offset_left
      const increment = p.seed_inc
      const increment_sm = 1
      curpos = p.seed_buffer
   else
      const increment = -p.seed_inc
      const increment_sm = -1
      curpos = readlen - p.seed_length - p.seed_buffer
   end
   const maxright = readlen - p.seed_length - p.seed_buffer
   while( ctry <= p.seed_try && p.seed_buffer <= curpos <= maxright )
      if minimum( read.metadata.quality[curpos:(curpos+p.seed_length-1)] ) <= p.seed_min_qual
         curpos += increment_sm
         continue
      end
      const curseq = read.seq[curpos:(curpos+p.seed_length-1)]
      const sa = FMIndexes.sa_range( curseq, index )
      const cnt = length(sa)
      ctry += 1
      #println(STDERR, "VERBOSE: SEED $sa, curpos: $curpos, cnt: $cnt, try: $(ctry-1)")
      if cnt == 0 || cnt > p.seed_tolerate
         if both_strands
            Bio.Seq.reverse_complement!(curseq)
            const rev_sa = FMIndexes.sa_range( curseq, index )
            const rev_cnt = length(rev_sa)
            if 0 < rev_cnt <= p.seed_tolerate
               const rev_pos = readlen - (curpos+p.seed_length-1)+1
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

function splice_by_identity!{A <: UngappedAlignment}( arr::Vector{A}, threshold, buffer, readlen )
   i = 1
   while i <= length( arr )
      if threshold - identity( arr[i], readlen ) > buffer
         splice!( arr, i )
         i -= 1
      end
      i += 1
   end
end

@inline function _ungapped_align( p::AlignParam, lib::GraphLib, read::SeqRecord, indx::Int, readloc::Int;
                                  ispos::Bool=true, geneind::NodeInt=convert(NodeInt,searchsortedlast( lib.offset, indx )) )

   align = ungapped_fwd_extend( p, lib, geneind,
                                indx - lib.offset[geneind] + 1,
                                read, readloc + 1, ispos=ispos )

   align = ungapped_rev_extend( p, lib, geneind,
                                indx - lib.offset[geneind],
                                read, readloc, ispos=ispos, align=align,
                                nodeidx=align.path[1].node )
   align
end

function ungapped_align( p::AlignParam, lib::GraphLib, read::SeqRecord; ispos::Bool=true, anchor_left::Bool=true )
   const seed,readloc,pos = seed_locate( p, lib.index, read, offset_left=anchor_left, both_strands=!p.is_stranded )
   const readlen = convert(UInt16, length(read.seq))

   if !pos
      Bio.Seq.reverse_complement!( read.seq )
      reverse!( read.metadata.quality )
      ispos = false
   end

   res      = Nullable{Vector{SGAlignment}}()
   maxscore = zero(Float32)

   for sidx in 1:length(seed)
      const s = FMIndexes.sa_value( seed[sidx], lib.index ) + 1

      const align = _ungapped_align( p, lib, read, s, readloc, ispos=ispos )
      const scvar = identity(align, readlen)

      if align.isvalid
         if isnull( res )
            res = Nullable(SGAlignment[ align ])
            maxscore = scvar
         else
            # new best score
            if scvar > maxscore
               # better than threshold for all of the previous scores
               if scvar - maxscore > p.score_range
                  length( res.value ) >= 1 && empty!( res.value )
               else # keep at least one previous score
                  splice_by_identity!( res.value, scvar, p.score_range, readlen ) 
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
   res
end


# This is the main ungapped alignment extension function in the --> direction
# Returns: SGAlignment
function ungapped_fwd_extend( p::AlignParam, lib::GraphLib, geneind::NodeInt, sgidx::Int, 
                                read::SeqRecord, ridx::Int; ispos::Bool=true,
                                align::SGAlignment=SGAlignment(0x0000, zero(Float32), sgidx, ridx, SGNode[],ispos,false),
                                nodeidx::NodeInt=convert(NodeInt,searchsortedlast(lib.graphs[geneind].nodeoffset,CoordInt(sgidx))) )
   const sg       = lib.graphs[geneind]
   const readlen  = convert(UInt16, length(read.seq))

   passed_edges = Nullable{Vector{UInt32}}() # don't allocate array unless needed
   passed_extend = 0
   passed_mismat = 0.0

   push!( align.path, SGNode( geneind, nodeidx ) ) # starting node

   #println(STDERR, "VERBOSE: FWD EXTENSION")

   while( align.mismatches <= p.mismatches && ridx <= readlen && sgidx <= length(sg.seq) )
      if sgidx == sg.nodeoffset[nodeidx] + sg.nodelen[nodeidx] # hit next edge
         # DEPRECATED isambiguous(sg.seq[sgidx]) && !(sg.seq[sgidx] == DNA_N) # L,R,S
         const curedge = nodeidx + 0x01
         if     sg.edgetype[curedge] == EDGETYPE_LR &&
                sg.nodelen[curedge]  >= p.kmer_size && # 'LR' && nodelen >= K
                readlen - ridx + 1   >= p.kmer_size
               # check edgeright[curnode+1] == read[ nextkmer ]
               # move forward K, continue 
               # This is a std exon-exon junction, if the read is an inclusion read
               # then we simply jump ahead, otherwise lets try to find a compatible right
               # node
               #println(STDERR, "VERBOSE: LR[1]")
               const rkmer_ind = kmer_index(read.seq[ridx:(ridx+p.kmer_size-1)])
               rkmer_ind == 0 && break
               if kmer_index(sg.edgeright[curedge]) == rkmer_ind
                  align.matches += convert(UInt16, p.kmer_size)
                  ridx    += p.kmer_size - 1
                  sgidx   += p.kmer_size - 1
                  nodeidx += 0x01
                  push!( align.path, SGNode( geneind, nodeidx ) )
                  align.isvalid = true
               else
                  align = spliced_fwd_extend( p, lib, geneind, curedge, read, ridx, rkmer_ind, align )
                  break
               end
         elseif sg.edgetype[curedge] == EDGETYPE_LR || 
                sg.edgetype[curedge] == EDGETYPE_LL # 'LR' || 'LL'
               # if length(edges) > 0, push! edgematch then set to 0
               # push! edge,  (here we try to extend fwd, but keep track
               # of potential edges we pass along the way.  When the alignment
               # dies if we have not had sufficient matches we then explore
               # those edges. 
               #println(STDERR, "VERBOSE: LR || LL")
               if isnull(passed_edges) # now we have to malloc
                  passed_edges = Nullable(Vector{UInt32}())
               end
               passed_extend = 0
               passed_mismat = 0.0
               push!(get(passed_edges), curedge)
               #sgidx   += 1
               #ridx    -= 1
               nodeidx += 0x01
               push!( align.path, SGNode( geneind, nodeidx ) )
         elseif sg.edgetype[curedge] == EDGETYPE_LS && # 'LS'
                readlen - ridx + 1   >= p.kmer_size
               # obligate spliced_extension
               const rkmer_ind = kmer_index(read.seq[ridx:(ridx+p.kmer_size-1)])
               rkmer_ind == 0 && break
               align = spliced_fwd_extend( p, lib, geneind, curedge, read, ridx, rkmer_ind, align )
               break
         elseif sg.edgetype[curedge] == EDGETYPE_SR #||
                #@bp
               # end of alignment
               break # ?
         else #'RR' || 'SL' || 'RS'
               # ignore 'RR' and 'SL'
               #sgidx += 1
               #ridx  -= 1 # offset the lower ridx += 1
               nodeidx += 0x01
               nodeidx <= length(sg.nodelen) && push!( align.path, SGNode( geneind, nodeidx ) )  
         end
         # ignore 'N'
      end
      if read.seq[ridx] == sg.seq[sgidx]
         # match
         align.matches += 0x01
         passed_extend += 1
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

   println(passed_edges)   
   println(passed_extend)
   # if edgemat < K, spliced_extension for each in length(edges)
   if !isnull(passed_edges)
      if passed_extend <= p.kmer_size + p.mismatches
         # go back.
         println("true")
         const ext_len = UInt16[sg.nodelen[ i ] for i in get(passed_edges)]
         ext_len[end] = passed_extend
         rev_cumarray!(ext_len)
         for c in length(get(passed_edges)):-1:1  #most recent edge first
            #println(STDERR, "VERBOSE: PASSED_EDGES - $c")
            if ext_len[c] >= p.kmer_size + p.mismatches
               align.isvalid = true
               break
            else
               pop!( align.path ) # extension into current node failed
               @fastmath align.mismatches -= passed_mismat
               align.matches -= ext_len[c]
               const cur_ridx = ridx - (sgidx - sg.nodeoffset[ get(passed_edges)[c] ])
               (cur_ridx + p.kmer_size - 1) <= readlen || continue
            
               #lkmer = DNAKmer{p.kmer_size}(read.seq[(ridx-p.kmer_size):(ridx-1)])
               const rkmer_ind  = kmer_index(read.seq[cur_ridx:(cur_ridx+p.kmer_size-1)])
               rkmer_ind == 0 && continue
               res_align = spliced_fwd_extend( p, lib, geneind, get(passed_edges)[c], 
                                               read, cur_ridx, rkmer_ind, align )
               #println(STDERR, "VERBOSE: RECURSIVE ALIGNMENT RESULT $res_align > $align $(res_align > align)")
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
   if identity(align, readlen) >= p.score_min && length(align.path) == 1  
      align.isvalid = true 
   end
   #println(STDERR, "VERBOSE: RETURNING $align")
   align
end



@inline function spliced_fwd_extend( p::AlignParam, lib::GraphLib, geneind::CoordInt, edgeind::UInt32, 
                             read::SeqRecord, ridx::Int, right_kmer_ind::Int, align::SGAlignment )
   # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
   # returns right nodes with matching genes
   isdefined( lib.edges.right, right_kmer_ind ) || return align

   best = align
   const left_kmer_ind = kmer_index(lib.graphs[geneind].edgeleft[edgeind])
   isdefined( lib.edges.left, left_kmer_ind ) || return align
   const right_nodes = lib.edges.left[ left_kmer_ind ] ∩ lib.edges.right[ right_kmer_ind ]

   #println(STDERR, "VERBOSE: SPLICED_EXTENSION, $(kmer_index(lib.graphs[geneind].edgeleft[edgeind])) ∩ $(kmer_index(rkmer)) = $right_nodes")

   for rn in right_nodes
      rn.gene == align.path[1].gene || continue
      const rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node]

      res_align::SGAlignment = ungapped_fwd_extend( p, lib, rn.gene, Int(rn_offset), read, ridx, 
                                                    align=deepcopy(align), nodeidx=rn.node )
      println(res_align)
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
function ungapped_rev_extend( p::AlignParam, lib::GraphLib, geneind::NodeInt, sgidx::Int, 
                                read::SeqRecord, ridx::Int; ispos::Bool=true,
                                align::SGAlignment=SGAlignment(0x0000, zero(Float32), sgidx, ridx, SGNode[], ispos, false),
                                nodeidx::NodeInt=convert(NodeInt,searchsortedlast(lib.graphs[geneind].nodeoffset,CoordInt(sgidx))) )
   const sg       = lib.graphs[geneind]
   const readlen  = convert(UInt16,length(read.seq))

   passed_edges = Nullable{Vector{UInt32}}() # don't allocate array unless needed
   passed_extend = 0
   passed_mismat = 0.0

   if align.path[1] != SGNode( geneind, nodeidx )
      unshift!( align.path, SGNode( geneind, nodeidx ) ) # starting node if not already there
   end

   while( align.mismatches <= p.mismatches && ridx > 0 && sgidx > 0 )
      if sgidx == sg.nodeoffset[nodeidx] - 1 # hit right edge
        # DEPRECATED isambiguous(sg.seq[sgidx]) && !(sg.seq[sgidx] == DNA_N) # L,R,S
         const leftnode = nodeidx - 0x01
         if     sg.edgetype[nodeidx] == EDGETYPE_LR &&
                sg.nodelen[leftnode] >= p.kmer_size && # 'LR' && nodelen >= K
                                ridx >= p.kmer_size 
                
               const lkmer_ind  = kmer_index(read.seq[(ridx-p.kmer_size+1):ridx])
               lkmer_ind == 0 && break 
               if kmer_index(sg.edgeleft[nodeidx]) == lkmer_ind
                  align.matches += convert(UInt16, p.kmer_size)
                  ridx    -= p.kmer_size - 1 
                  sgidx   -= p.kmer_size - 1 
                  nodeidx -= 0x01
                  unshift!( align.path, SGNode( geneind, nodeidx ) )
                  align.isvalid = true
               else
                  align = spliced_rev_extend( p, lib, geneind, nodeidx, read, ridx, lkmer_ind, align ) 
                  break
               end
         elseif sg.edgetype[nodeidx] == EDGETYPE_LR || 
                sg.edgetype[nodeidx] == EDGETYPE_RR # 'LR' || 'RR'
                
               if isnull(passed_edges) # now we have to malloc
                  passed_edges = Nullable(Vector{UInt32}())
               end
               passed_extend = 0
               passed_mismat = 0.0
               push!(get(passed_edges), nodeidx)
               #sgidx   -= 1
               #ridx    += 1
               nodeidx -= 0x01
               unshift!( align.path, SGNode( geneind, nodeidx ) )
         elseif sg.edgetype[nodeidx] == EDGETYPE_SR && # 'SR'
                                ridx >= p.kmer_size
               # obligate spliced_extension
               const lkmer_ind  = kmer_index(read.seq[(ridx-p.kmer_size+1):ridx])
               lkmer_ind == 0 && break
               align = spliced_rev_extend( p, lib, geneind, nodeidx, read, ridx, lkmer_ind, align )
               break
         elseif sg.edgetype[nodeidx] == EDGETYPE_LS #||
               # end of alignment
               break # ?
         else #'LL' || 'RS' || 'SL'
               # ignore 'LL' and 'RS'
               #sgidx -= 1
               #ridx  += 1 # offset the lower ridx += 1
               nodeidx -= 0x01
               nodeidx > 0 && unshift!( align.path, SGNode( geneind, nodeidx ) )  
         end
         # ignore 'N'
      end
      if read.seq[ridx] == sg.seq[sgidx]
         # match
         align.matches += 0x01
         passed_extend += 1 
         # mismatch
      else
         const prob  = phred_to_prob( read.metadata.quality[ridx] )
         @fastmath align.mismatches += prob
         @fastmath passed_mismat += prob
      end
      ridx  -= 1
      sgidx -= 1
      #print(STDERR, " $(read.seq[ridx+1]),$ridx\_$(sg.seq[sgidx+1]),$sgidx ")
   end

   cur_ridx  = ridx + 1
   cur_sgidx = sgidx + 1
   # if passed_extend < K, spliced_extension for each in length(edges)
   if !isnull(passed_edges)
      if passed_extend <= p.kmer_size + p.mismatches
         # go back.
         const ext_len = UInt16[sg.nodelen[ i ] for i in get(passed_edges)]
         ext_len[end] = passed_extend
         rev_cumarray!(ext_len)
         for c in length(get(passed_edges)):-1:1  #most recent edge first
            if ext_len[c] >= p.kmer_size
               align.isvalid = true
               break
            else
               shift!( align.path ) # extension into current node failed
               @fastmath align.mismatches -= passed_mismat
               align.matches -= ext_len[c]
               cur_ridx  = ridx + (sg.nodeoffset[ get(passed_edges)[c] ] - sgidx)
               cur_sgidx = sg.nodeoffset[ get(passed_edges)[c] ]
               (cur_ridx - p.kmer_size) > 0 || continue
            
               const lkmer_ind  = kmer_index(read.seq[(cur_ridx-p.kmer_size):(cur_ridx-1)])
               lkmer_ind == 0 && continue
               res_align = spliced_rev_extend( p, lib, geneind, get(passed_edges)[c], 
                                                  read, cur_ridx-1, lkmer_ind, align )
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
   if identity(align, readlen) >= p.score_min && length(align.path) == 1  
      align.isvalid = true 
   end
   if sgidx < align.offset
      align.offset = convert(UInt32, cur_sgidx)
   end
   if cur_ridx < align.offsetread
      align.offsetread = convert(UInt16, cur_ridx)
   end
   align
end

@inline function spliced_rev_extend( p::AlignParam, lib::GraphLib, geneind::NodeInt, edgeind::CoordInt,
                             read::SeqRecord, ridx::Int, left_kmer_index::Int, align::SGAlignment )
   # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
   # returns right nodes with matching genes
   isdefined( lib.edges.left, left_kmer_index ) || return align

   best = align
   const right_kmer_ind = kmer_index(lib.graphs[geneind].edgeright[edgeind])
   isdefined( lib.edges.right, right_kmer_ind ) || return align
   const left_nodes = lib.edges.right[ right_kmer_ind ] ∩ lib.edges.left[ left_kmer_index ]

   # do a test for trans-splicing, and reset left_nodes
   for rn in left_nodes
      rn.gene == align.path[1].gene || continue
      const rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node-1] + lib.graphs[rn.gene].nodelen[rn.node-1] - 1

      res_align::SGAlignment = ungapped_rev_extend( p, lib, rn.gene, Int(rn_offset), read, ridx, 
                                                    align=deepcopy(align), nodeidx=rn.node-0x01 )
      best = res_align > best ? res_align : best
   end
   if score(best) >= score(align) + p.kmer_size
      best.isvalid = true
   end
   best
end
