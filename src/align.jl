
struct AlignParam
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
@inline function AlignParam( args::Dict{S,Any}, ispaired=false; kmer=9 ) where S <: AbstractString
   aln = AlignParam( args["mismatches"], kmer,
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

abstract type UngappedAlignment end

const MatchesInt    = UInt8
const MismatchInt   = UInt8
const MismatchFloat = Float32

mutable struct SGAlignScore
   matches::MatchesInt
   mismatches::MismatchInt
   mistolerance::MismatchFloat
end

@inline Base.zero(t::Type{SGAlignScore}) = SGAlignScore(0x0000,zero(MismatchInt),zero(MismatchFloat))
@inline Base.one(t::Type{SGAlignScore}) = SGAlignScore(0x0001,one(MismatchInt),one(MismatchFloat))
@inline match(t::Type{SGAlignScore}, val::I) where {I <: Integer} = SGAlignScore(convert(MatchesInt, val),zero(MismatchInt),zero(MismatchFloat))

struct SGAlignNode
   gene::NodeInt
   node::NodeInt
   score::SGAlignScore
end

Base.hash( sgn::SGAlignNode ) = hash( sgn.gene, hash(sgn.node) )
Base.hash( sgn::SGAlignNode, h::UInt64 ) = hash( sgn.gene, hash(sgn.node, h) )
Base.isequal( a::SGAlignNode, b::SGAlignNode ) = a.gene == b.gene && a.node == b.node
Base.:(==)( a::SGAlignNode, b::SGAlignNode ) = a.gene == b.gene && a.node == b.node

const SGAlignPath = Vector{SGAlignNode}

@inline SGAlignNode(sgn::SGNode) = SGAlignNode(sgn.gene, sgn.node, zero(SGAlignScore))
@inline Base.push!( v::Vector{SGAlignNode}, sgn::SGNode ) = push!( v, SGAlignNode(sgn) )

const ReadLengthInt = UInt8

mutable struct SGAlignment <: UngappedAlignment
   offset::UInt32
   offsetread::ReadLengthInt
   path::Vector{SGAlignNode}
   strand::Bool
   isvalid::Bool
end

SGAlignment() = SGAlignment(zero(UInt32), zero(ReadLengthInt), SGAlignNode[], true, false)

const SGAlignVec = Nullable{Vector{SGAlignment}}

const DEF_ALIGN = SGAlignment(zero(UInt32), zero(ReadLengthInt), SGAlignNode[], true, false)

@inbounds function matches( v::Vector{SGAlignNode} )
   matches = zero(MatchesInt)
   for i in 1:length(v)
      matches += v[i].score.matches
   end
   matches
end
@inbounds function mistolerance( v::Vector{SGAlignNode} )
   mistol = zero(MismatchFloat)
   for i in 1:length(v)
      mistol += v[i].score.mistolerance
   end
   mistol
end
@inbounds function score( v::Vector{SGAlignNode} )
   score = zero(MismatchFloat)
   for i in 1:length(v)
      score += v[i].score.matches - v[i].score.mistolerance
   end
   score
end
@inline matches( aln::A ) where {A <: UngappedAlignment} = matches( aln.path )
@inline mistolerance( aln::A ) where {A <: UngappedAlignment} = mistolerance( aln.path )
@inline score( aln::A ) where {A <: UngappedAlignment} = score( aln.path )
@inline identity( align::A, readlen::I ) where {A <: UngappedAlignment, I <: Integer} = convert(Float32, @fastmath score( align ) / readlen)
@inline isvalid( align::A ) where {A <: UngappedAlignment} = align.isvalid && length(align.path) >= 1 && matches(align) > mistolerance(align) ? true : false

Base.:>( a::SGAlignment, b::SGAlignment ) = >( score(a), score(b) )
Base.:<( a::SGAlignment, b::SGAlignment ) = <( score(a), score(b) )
Base.:>=( a::SGAlignment, b::SGAlignment ) = >=( score(a), score(b) )
Base.:<=( a::SGAlignment, b::SGAlignment ) = <=( score(a), score(b) )
Base.isless( a::SGAlignment, b::SGAlignment ) = a < b

# add prob of being accurate base to mismatch, rather than integer.
@inline phred_to_prob( phred::UInt8 ) = convert(Float32, @fastmath 1-10^(-convert(Int8,phred)/10))

@inline function Base.empty!( align::SGAlignment )
   align.offset = zero(UInt32)
   align.offsetread = zero(ReadLengthInt)
   empty!( align.path )
   align.strand = true
   align.isvalid = false
end

@inline function remove_invalid!( vec::Vector{SGAlignment} )
   i = 1
   while i <= length(vec)
      @inbounds if !isvalid( vec[i] )
         deleteat!( vec, i )
      else
         i += 1
      end
   end
end

@inline function seed_locate( p::AlignParam, index::FMIndex, read::FASTQRecord; offset_left::Bool=true, both_strands::Bool=true )
   def_sa = 2:1
   readlen = convert(ReadLengthInt, length(read.sequence))
   ctry = 1
   if offset_left
      increment = p.seed_inc
      increment_sm = 1
      curpos = p.seed_buffer
   else
      increment = -p.seed_inc
      increment_sm = -1
      curpos = readlen - p.seed_length - p.seed_buffer
   end
   maxright = readlen - p.seed_length - p.seed_buffer
   while( ctry <= p.seed_try && p.seed_buffer <= curpos <= maxright )
      if minimum( read.quality[curpos:(curpos+p.seed_length-1)] ) <= p.seed_min_qual
         curpos += increment_sm
         continue
      end
      curseq = read.sequence[curpos:(curpos+p.seed_length-1)]
      sa = FMIndexes.sa_range( curseq, index )
      cnt = length(sa)
      ctry += 1
      if cnt == 0 || cnt > p.seed_tolerate
         if both_strands
            BioSequences.reverse_complement!(curseq)
            # curseq = reverse_complement(curseq)
            rev_sa = FMIndexes.sa_range( curseq, index )
            rev_cnt = length(rev_sa)
            if 0 < rev_cnt <= p.seed_tolerate
               rev_pos = readlen - (curpos+p.seed_length-1)+1
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

@inline function splice_by_score!( arr::Vector{A}, threshold, buffer ) where A <: UngappedAlignment
   i = 1
   while i <= length( arr )
      @inbounds if threshold - score( arr[i] ) > buffer
         splice!( arr, i )
         i -= 1
      end
      i += 1
   end
end

@inline function splice_by_identity!( arr::Vector{A}, threshold, buffer, readlen ) where A <: UngappedAlignment
   i = 1
   while i <= length( arr )
      @inbounds if threshold - identity( arr[i], readlen ) > buffer
         splice!( arr, i )
         i -= 1
      end
      i += 1
   end
end

@inline function _ungapped_align( p::AlignParam, lib::GraphLib, read::FASTQRecord, indx::Int, readloc::Int;
                                  ispos::Bool=true, geneind::NodeInt=convert(NodeInt,searchsortedlast( lib.offset, indx )) )

   sgidx       = indx - lib.offset[geneind]
   offset_node = convert(NodeInt,searchsortedlast(lib.graphs[geneind].nodeoffset,convert(CoordInt, sgidx + 1)))

   align = ungapped_fwd_extend( p, lib, geneind, sgidx + 1,
                                read, readloc + 1, ispos=ispos,
                                nodeidx=offset_node )

   align = ungapped_rev_extend( p, lib, geneind, sgidx,
                                read, readloc, ispos=ispos, align=align,
                                nodeidx=offset_node )
   align
end

@inline function ungapped_align( p::AlignParam, lib::GraphLib, read::FASTQRecord;
                         ispos::Bool=true, anchor_left::Bool=true )

   seed,readloc,pos = seed_locate( p, lib.index, read, offset_left=anchor_left, both_strands=!p.is_stranded )
   readlen = convert(ReadLengthInt, length(read.sequence))

   if !pos
      reverse_complement!( read )
      ispos = false
   end

   res      = Nullable{Vector{SGAlignment}}()
   maxscore = zero(Float32)

   for sidx in 1:length(seed)
      s = FMIndexes.sa_value( seed[sidx], lib.index ) + 1

      align = _ungapped_align( p, lib, read, s, readloc, ispos=ispos )
      scvar = identity(align, readlen)

      if isvalid(align)
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
@inbounds function ungapped_fwd_extend( p::AlignParam, lib::GraphLib, geneind::NodeInt, sgidx::Int,
                                read::FASTQRecord, ridx::Int; ispos::Bool=true,
                                align::SGAlignment=SGAlignment(sgidx, ridx, SGAlignNode[],ispos,false),
                                nodeidx::NodeInt=convert(NodeInt,searchsortedlast(lib.graphs[geneind].nodeoffset,CoordInt(sgidx))) )
   sg       = lib.graphs[geneind]
   readlen  = convert(ReadLengthInt, length(read.sequence))

   passed_edges  = Nullable{Vector{UInt8}}() # don't allocate array unless needed
   mistol        = mistolerance(align)

   push!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) ) # starting node

   @inline function spliced_fwd_extend( geneind::CoordInt, edgeind::UInt32,
                                        ridx::Int, right_kmer_ind::Int, align::SGAlignment )
      # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
      # returns right nodes with matching genes
      isassigned( lib.edges.right, right_kmer_ind ) || return align

      best = align
      left_kmer_ind = kmer_index(lib.graphs[geneind].edgeleft[edgeind])
      isassigned( lib.edges.left, left_kmer_ind ) || return align
      @inbounds right_nodes = lib.edges.left[ left_kmer_ind ] ∩ lib.edges.right[ right_kmer_ind ]

      for rn in right_nodes
         rn.gene == align.path[1].gene || continue
         rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node]

         res_align::SGAlignment = ungapped_fwd_extend( p, lib, rn.gene, Int(rn_offset), read, ridx,
                                                       align=deepcopy(align), nodeidx=rn.node )
         #println(res_align)
         best = res_align > best ? res_align : best
      end
      if score(best) >= score(align) + p.kmer_size
         best.isvalid = true
      end
      best
   end

   while( mistol <= p.mismatches && ridx <= readlen && sgidx <= length(sg.seq) )
      if sgidx == sg.nodeoffset[nodeidx] + sg.nodelen[nodeidx] # hit next edge
         curedge = nodeidx + 0x01
         if     sg.edgetype[curedge] == EDGETYPE_LR &&
                sg.nodelen[curedge]  >= p.kmer_size && # 'LR' && nodelen >= K
                readlen - ridx + 1   >= p.kmer_size
               # check edgeright[curnode+1] == read[ nextkmer ]
               # move forward K, continue
               # This is a std exon-exon junction, if the read is an inclusion read
               # then we simply jump ahead, otherwise lets try to find a compatible right
               # node
               rkmer_ind = kmer_index(read.sequence[ridx:(ridx+p.kmer_size-1)])
               rkmer_ind == 0 && break
               if kmer_index(sg.edgeright[curedge]) == rkmer_ind
                  ridx    += p.kmer_size - 1
                  sgidx   += p.kmer_size - 1
                  nodeidx += 0x01
                  push!( align.path, SGAlignNode( geneind, nodeidx, match(SGAlignScore, p.kmer_size - 1) ) )
                  align.isvalid = true
               else
                  align = spliced_fwd_extend( geneind, curedge, ridx, rkmer_ind, align )
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
                  passed_edges = Nullable(Vector{UInt8}())
               end
               push!(get(passed_edges), length(align.path))
               nodeidx += 0x01
               push!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) )
         elseif sg.edgetype[curedge] == EDGETYPE_LS # 'LS'
               if readlen - ridx + 1   >= p.kmer_size
                  # obligate spliced_extension
                  rkmer_ind = kmer_index(read.sequence[ridx:(ridx+p.kmer_size-1)])
                  rkmer_ind == 0 && break
                  align = spliced_fwd_extend( geneind, curedge, ridx, rkmer_ind, align )
               end
               break
         elseif sg.edgetype[curedge] == EDGETYPE_SR
               # end of alignment
               break # ?
         else #'RR' || 'SL' || 'RS'
               # ignore 'RR' and 'SL'
               nodeidx += 0x01
               nodeidx <= length(sg.nodelen) && push!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) )
         end
         # ignore 'N'
      end
      if read.sequence[ridx] == sg.seq[sgidx]
         # match
         align.path[end].score.matches += 0x01
      else
         # mismatch
         prob  = phred_to_prob( read.quality[ridx] )
         @fastmath align.path[end].score.mistolerance += prob
         @fastmath mistol += prob
         align.path[end].score.mismatches += 0x01
      end
      ridx  += 1
      sgidx += 1
   end

   # if edgemat < K, spliced_extension for each in length(edges)
   if !isnull(passed_edges)
         # go back.
      for c in length(get(passed_edges)):-1:1  #most recent edge first
         cval = get(passed_edges)[c]
         if matches(align.path[cval+1:end]) >= p.kmer_size + p.mismatches
            align.isvalid = true
            break
         else
            cur_edge = align.path[cval].node + one(UInt32)
            splice!( align.path, (cval+1):length(align.path) ) # extension into current node failed

            cur_ridx = ridx - (sgidx - sg.nodeoffset[ cur_edge ])
            (cur_ridx + p.kmer_size - 1) <= readlen || continue

            rkmer_ind  = kmer_index(read.sequence[cur_ridx:(cur_ridx+p.kmer_size-1)])
            rkmer_ind == 0 && continue
            res_align = spliced_fwd_extend( geneind, cur_edge,
                                            cur_ridx, rkmer_ind, align )
            if length(res_align.path) > length(align.path) # we added a valid node
               align = res_align > align ? res_align : align
            end
         end
      end
   end

   # clean up any bad trailing nodes
   while length(align.path) >= 1 && (align.path[end].score.mismatches >= align.path[end].score.matches ||
                                     Int(align.path[end].score.matches) - Int(align.path[end].score.mismatches) <= 1)
      pop!( align.path )
      length(align.path) == 0 && (align.isvalid = false)
   end

   if identity(align, readlen) >= p.score_min
      align.isvalid = true
   end
   align
end


# This is the main ungapped alignment extension function in the <-- direction
# Returns: SGAlignment
@inbounds function ungapped_rev_extend( p::AlignParam, lib::GraphLib, geneind::NodeInt, sgidx::Int,
                                read::FASTQRecord, ridx::Int; ispos::Bool=true,
                                align::SGAlignment=SGAlignment(sgidx, ridx, SGAlignNode[], ispos, false),
                                nodeidx::NodeInt=convert(NodeInt,searchsortedlast(lib.graphs[geneind].nodeoffset,CoordInt(sgidx))) )
   sg       = lib.graphs[geneind]
   readlen  = convert(ReadLengthInt,length(read.sequence))

   passed_edges  = Nullable{Vector{UInt8}}() # don't allocate array unless needed
   mistol        = mistolerance(align)

   if length(align.path) <= 0 || align.path[1].gene != geneind || align.path[1].node != nodeidx
      pushfirst!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) ) # starting node if not already there
   end

   @inline function spliced_rev_extend( geneind::NodeInt, edgeind::CoordInt,
                                        ridx::Int, left_kmer_index::Int, align::SGAlignment )
      # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
      # returns right nodes with matching genes
      isassigned( lib.edges.left, left_kmer_index ) || return align

      best = align
      right_kmer_ind = kmer_index(lib.graphs[geneind].edgeright[edgeind])
      isassigned( lib.edges.right, right_kmer_ind ) || return align
      left_nodes = lib.edges.right[ right_kmer_ind ] ∩ lib.edges.left[ left_kmer_index ]

      # do a test for trans-splicing, and reset left_nodes
      for rn in left_nodes
         rn.gene == align.path[1].gene || continue
         rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node-1] + lib.graphs[rn.gene].nodelen[rn.node-1] - 1

         res_align::SGAlignment = ungapped_rev_extend( p, lib, rn.gene, Int(rn_offset), read, ridx,
                                                       align=deepcopy(align), nodeidx=rn.node-0x01 )
         best = res_align > best ? res_align : best
      end
      if score(best) >= score(align) + p.kmer_size
         best.isvalid = true
      end
      best
   end

   while( mistol <= p.mismatches && ridx > 0 && sgidx > 0 )
      if sgidx == sg.nodeoffset[nodeidx] - 1 # hit right edge
         leftnode = nodeidx - 0x01
         if     sg.edgetype[nodeidx] == EDGETYPE_LR &&
                sg.nodelen[leftnode] >= p.kmer_size && # 'LR' && nodelen >= K
                                ridx >= p.kmer_size

               lkmer_ind  = kmer_index(read.sequence[(ridx-p.kmer_size+1):ridx])
               lkmer_ind == 0 && break
               if kmer_index(sg.edgeleft[nodeidx]) == lkmer_ind
                  ridx    -= p.kmer_size - 1
                  sgidx   -= p.kmer_size - 1
                  nodeidx -= 0x01
                  pushfirst!( align.path, SGAlignNode( geneind, nodeidx, match(SGAlignScore, p.kmer_size - 1) ) )
                  align.isvalid = true
               else
                  align = spliced_rev_extend( geneind, nodeidx, ridx, lkmer_ind, align )
                  break
               end
         elseif sg.edgetype[nodeidx] == EDGETYPE_LR ||
                sg.edgetype[nodeidx] == EDGETYPE_RR # 'LR' || 'RR'

               if isnull(passed_edges) # now we have to malloc
                  passed_edges = Nullable(Vector{UInt8}())
               end
               push!(get(passed_edges), length(align.path)-1)
               nodeidx -= 0x01
               pushfirst!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) )
         elseif sg.edgetype[nodeidx] == EDGETYPE_SR # 'SR'
               if ridx >= p.kmer_size
                  # obligate spliced_extension
                  lkmer_ind  = kmer_index(read.sequence[(ridx-p.kmer_size+1):ridx])
                  lkmer_ind == 0 && break
                  align = spliced_rev_extend( geneind, nodeidx, ridx, lkmer_ind, align )
               end
               break
         elseif sg.edgetype[nodeidx] == EDGETYPE_LS
               # end of alignment
               break # ?
         else #'LL' || 'RS' || 'SL'
               # ignore 'LL' and 'RS'
               nodeidx -= 0x01
               nodeidx > 0 && pushfirst!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) )
         end
         # ignore 'N'
      end
      if read.sequence[ridx] == sg.seq[sgidx]
         # match
         align.path[1].score.matches += 0x01
         # mismatch
      else
         prob  = phred_to_prob( read.quality[ridx] )
         @fastmath align.path[1].score.mistolerance += prob
         @fastmath mistol += prob
         align.path[1].score.mismatches += 0x01
      end
      ridx  -= 1
      sgidx -= 1
   end

   cur_ridx  = ridx + 1
   cur_sgidx = sgidx + 1
   # if passed_extend < K, spliced_extension for each in length(edges)
   if !isnull(passed_edges)
         # go back.
      for c in length(get(passed_edges)):-1:1  #most recent edge first
         cval = get(passed_edges)[c]
         if matches(align.path[1:(end-(cval+1))]) >= p.kmer_size + p.mismatches
            align.isvalid = true
            break
         else
            cur_edge = align.path[end-cval].node
            splice!( align.path, 1:(length(align.path)-(cval+1)) ) # extension into current node failed
            cur_ridx  = ridx + (sg.nodeoffset[ cur_edge ] - sgidx)
            cur_sgidx = sg.nodeoffset[ cur_edge ]
            (cur_ridx - p.kmer_size) > 0 || continue

            lkmer_ind  = kmer_index(read.sequence[(cur_ridx-p.kmer_size):(cur_ridx-1)])
            lkmer_ind == 0 && continue
            res_align = spliced_rev_extend( geneind, cur_edge,
                                            cur_ridx-1, lkmer_ind, align )
            if length(res_align.path) > length(align.path) # we added a valid node
               align = res_align > align ? res_align : align
            end
         end
      end
   end

   if identity(align, readlen) >= p.score_min
      align.isvalid = true
   end
   if sgidx < align.offset
      align.offset = convert(UInt32, cur_sgidx)
   end
   if cur_ridx < align.offsetread
      align.offsetread = convert(ReadLengthInt, cur_ridx)
   end

   # clean up any bad left trailing nodes
   while length(align.path) >= 1 && (align.path[1].score.mismatches >= align.path[1].score.matches ||
                                     Int(align.path[1].score.matches) - Int(align.path[1].score.mismatches) <= 1)
      offset_change = align.path[1].score.matches + align.path[1].score.mismatches
      popfirst!( align.path )
      if length(align.path) >= 1
         align.offset      = sg.nodeoffset[ align.path[1].node ]
         align.offsetread += offset_change
      else
         align.isvalid = false
      end
   end

   align
end
