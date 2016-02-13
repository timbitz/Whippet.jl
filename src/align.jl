# requires
# using FMIndexes
# include("graph.jl")
# include("index.jl")

import Base.<,Base.>,Base.<=,Base.>=

import Bio.Seq.SeqRecord

immutable AlignParam
   mismatches::Int    # Allowable mismatches
   kmer_size::Int     # Minimum number of matches to extend past an edge
   seed_try::Int      # Starting number of seeds
   seed_tolerate::Int # Allow at most _ hits for a valid seed
   seed_length::Int   # Seed Size
   seed_maxoff::Int   # Seed gathering constraint
   seed_buffer::Int   # Ignore first _ bases
   seed_inc::Int      # Incrementation for subsequent seed searches
   seed_rng::Int      # Seed for random number generator
   score_range::Int   # Scoring range to be considered repetitive
   score_min::Int     # Minimum score for a valid alignment
   is_stranded::Bool  # Is input data + strand only?
   is_paired::Bool    # Paired end data?
   is_trans_ok::Bool  # Do we allow edges from one gene to another
   is_circ_ok::Bool   # Do we allow back edges
end

AlignParam() = AlignParam( 2, 9, 4, 7, 18, 35, 5, 10, 1, 10, 25, false, false, false, true )

abstract UngappedAlignment

type SGAlignment <: UngappedAlignment
   matches::Int
   mismatches::Float64
   offset::Int
   path::Vector{SGNode}
   strand::Bool
   isvalid::Bool
end

const DEF_ALIGN = SGAlignment(0, 0, 0, SGNode[], true, false)

score( align::SGAlignment ) = align.matches - align.mismatches 

>( a::SGAlignment, b::SGAlignment ) = >( score(a), score(b) )
<( a::SGAlignment, b::SGAlignment ) = <( score(a), score(b) )
>=( a::SGAlignment, b::SGAlignment ) = >=( score(a), score(b) )
<=( a::SGAlignment, b::SGAlignment ) = <=( score(a), score(b) )

# add prob of being accurate base to mismatch, rather than integer.
phred_to_prob( phred::Int8 ) = @fastmath 1-10^(-phred/10)

function seed_locate( p::AlignParam, index::FMIndex, read::SeqRecord )
   sa = 2:1
   cnt = 1000
   ctry = 0
   curpos = p.seed_buffer - p.seed_inc
   while( (cnt == 0 || cnt > p.seed_tolerate) && 
          ctry <= p.seed_try && curpos <= p.seed_maxoff )
      curpos += p.seed_inc
      sa = FMIndexes.sa_range( read.seq[curpos:(curpos+p.seed_length-1)], index )
      cnt = length(sa)
      ctry += 1
   end
   sa,curpos
end

function ungapped_align( p::AlignParam, lib::GraphLib, read::SeqRecord; ispos=true )
   seed,readloc = seed_locate( p, lib.index, read )
   locit = FMIndexes.LocationIterator( seed, lib.index )
   #@bp
   res   = Nullable{Vector{SGAlignment}}()
   for s in locit
      geneind = search_sorted( lib.offset, convert(Coordint, s), lower=true ) 
      #println("$(read.seq[readloc:(readloc+75)])\n$(lib.graphs[geneind].seq[(s-lib.offset[geneind]):(s-lib.offset[geneind])+50])")
      #@bp
      align = ungapped_fwd_extend( p, lib, convert(Coordint, geneind), s - lib.offset[geneind] + p.seed_length, 
                                   read, readloc + p.seed_length, ispos=ispos ) # TODO check
 
      #align = ungapped_rev_extend( p, lib, s - lib.offset[geneid] - 1,
      #                             read, readloc - 1, ispos=ispos, align=align, nodeidx=align.path[1][2] )
      if align.isvalid
         if isnull( res )
            res = Nullable(Vector{SGAlignment}())
         end
         push!(get(res), align)
      end
   end
   # if !stranded and no valid alignments, run reverse complement
   if ispos && !p.is_stranded && isnull( res )
      read.seq = Bio.Seq.reverse_complement( read.seq )
      res = ungapped_align( p, lib, read, ispos=false )
   end
   res
end


# This is the main ungapped alignment extension function in the --> direction
# Returns: SGAlignment
function ungapped_fwd_extend( p::AlignParam, lib::GraphLib, geneind::Coordint, sgidx::Int, 
                              read::SeqRecord, ridx::Int; ispos=true,
                              align=SGAlignment(p.seed_length,0,sgidx,SGNode[],ispos,false),
                              nodeidx=search_sorted(lib.graphs[geneind].nodeoffset,Coordint(sgidx),lower=true) )
   sg       = lib.graphs[geneind]
   readlen  = length(read.seq)
   curedge  = nodeidx

   passed_edges = Nullable{Vector{Coordint}}() # don't allocate array unless needed
   passed_extend = 0
   passed_mismat = 0.0

   push!( align.path, SGNode( geneind, nodeidx ) ) # starting node

   while( align.mismatches <= p.mismatches && ridx <= readlen && sgidx <= length(sg.seq) )
      #@bp
      if read.seq[ridx] == sg.seq[sgidx]
         # match
         align.matches += 1
         passed_extend  += 1
      elseif (UInt8(sg.seq[sgidx]) & 0b100) == 0b100 && !(sg.seq[sgidx] == SG_N) # L,R,S
         curedge = nodeidx+1
         if     sg.edgetype[curedge] == EDGETYPE_LR &&
                sg.nodelen[curedge]  >= p.kmer_size && # 'LR' && nodelen >= K
                readlen - ridx + 1   >= p.kmer_size
               # check edgeright[curnode+1] == read[ nextkmer ]
               # move forward K, continue 
               # This is a std exon-exon junction, if the read is an inclusion read
               # then we simply jump ahead, otherwise lets try to find a compatible right
               # node
               rkmer = DNAKmer{p.kmer_size}(read.seq[ridx:(ridx+p.kmer_size-1)])
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
                  passed_edges = Nullable(Vector{Coordint}())
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
               rkmer = DNAKmer{p.kmer_size}(read.seq[ridx:(ridx+p.kmer_size-1)])
               align = spliced_fwd_extend( p, lib, geneind, curedge, read, ridx, rkmer, align )
               break
         elseif sg.edgetype[curedge] == EDGETYPE_SR ||
                sg.edgetype[curedge] == EDGETYPE_RS # 'SR' || 'RS'
               # end of alignment
               break # ?
         else #'RR' || 'SL'
               # ignore 'RR' and 'SL'
               sgidx += 1
               ridx  -= 1 # offset the lower ridx += 1
               nodeidx += 1
               push!( align.path, SGNode( geneind, nodeidx ) )  
         end
         # ignore 'N'
      else 
         # mismatch
         align.mismatches += phred_to_prob( read.metadata.quality[ridx] )
         passed_mismat += phred_to_prob( read.metadata.quality[ridx] )
      end
      ridx  += 1
      sgidx += 1
      #print(" $(read.seq[ridx-1]),$ridx\_$(sg.seq[sgidx-1]),$sgidx ")
   end

   # if edgemat < K, spliced_extension for each in length(edges)
   # TODO go through passed_edges!!
   if !isnull(passed_edges)
      if passed_extend < p.kmer_size
         # go back.
         ext_len = Int[sg.nodelen[ i ] for i in get(passed_edges)]
         ext_len[end] = passed_extend
         rev_cumarray!(ext_len)
         #@bp
         for c in length(get(passed_edges)):-1:1  #most recent edge first
            if ext_len[c] >= p.kmer_size
               align.isvalid = true
               break
            else
               pop!( align.path ) # extension into current node failed
               align.mismatches -= passed_mismat
               cur_ridx = ridx - (sgidx - sg.nodeoffset[ get(passed_edges)[c] ])
               (cur_ridx + p.kmer_size - 1) < readlen || continue
            
               #lkmer = DNAKmer{p.kmer_size}(read.seq[(ridx-p.kmer_size):(ridx-1)])
               rkmer = DNAKmer{p.kmer_size}(read.seq[cur_ridx:(cur_ridx+p.kmer_size-1)])
               #println("$(read.seq[(ridx-p.kmer_size):(ridx-1)])\n$lkmer\n$(sg.edgeleft[get(passed_edges)[c]])")

               res_align = spliced_fwd_extend( p, lib, geneind, get(passed_edges)[c], read, cur_ridx, rkmer, align )
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



function spliced_fwd_extend( p::AlignParam, lib::GraphLib, geneind, edgeind, read, ridx, rkmer, align )
   # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
   # returns right nodes with matching genes
   isdefined( lib.edges.right, kmer_index(rkmer) ) || return DEF_ALIGN

   best = DEF_ALIGN
   left_kmer = lib.graphs[geneind].edgeleft[edgeind]
   right_nodes = lib.edges.left[ kmer_index(left_kmer) ] ∩ lib.edges.right[ kmer_index(rkmer) ]

   # do a test for trans splicing, and reset right_nodes
   for rn in right_nodes
      rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node]
      res_align = ungapped_fwd_extend( p, lib, rn.gene, Int(rn_offset), read, ridx, 
                                       align=deepcopy(align), nodeidx=rn.node )
      best = res_align > best ? res_align : best
   end
   if score(best) >= score(align) + p.kmer_size
      best.isvalid = true
   end
   best
end

# Function Cumulative Array
function rev_cumarray!{T<:Vector}( array::T )
   for i in (length(array)-1):-1:1
      array[i] += array[i+1]
   end 
end


# This is the main ungapped alignment extension function in the <-- direction
# Returns: SGAlignment
function ungapped_rev_extend( p::AlignParam, lib::GraphLib, geneind::Coordint, sgidx::Int, 
                              read::SeqRecord, ridx::Int; ispos=true,
                              align=SGAlignment(p.seed_length,0,sgidx,SGNode[],ispos,false),
                              nodeidx=search_sorted(lib.graphs[geneind].nodeoffset,Coordint(sgidx),lower=true) )
   sg       = lib.graphs[geneind]
   readlen  = length(read.seq)
   curedge  = nodeidx

   passed_edges = Nullable{Vector{Coordint}}() # don't allocate array unless needed
   passed_extend = 0
   passed_mismat = 0

   if align.path[1] != SGNode( geneind, nodeidx )
      unshift!( align.path, SGNode( geneind, nodeidx ) ) # starting node if not already there
   end

   while( align.mismatches <= p.mismatches && ridx > 0 && sgidx > 0 )
      #@bp
      if read.seq[ridx] == sg.seq[sgidx]
         # match
         align.matches -= 1
         passed_extend  -= 1
      elseif (UInt8(sg.seq[sgidx]) & 0b100) == 0b100 && !(sg.seq[sgidx] == SG_N) # L,R,S
         leftnode = nodeidx - 1
         if     sg.edgetype[nodeidx] == EDGETYPE_LR &&
                sg.nodelen[leftnode] >= p.kmer_size && # 'LR' && nodelen >= K
                                ridx >= p.kmer_size  # TODO check
                                
               lkmer = DNAKmer{p.kmer_size}(read.seq[(ridx-p.kmer_size+1):ridx])
               if sg.edgeleft[nodeidx] == lkmer
                  align.matches += p.kmer_size
                  ridx    -= p.kmer_size + 1
                  sgidx   -= 1 + p.kmer_size # TODO check
                  nodeidx -= 1
                  unshift!( align.path, SGNode( geneind, nodeidx ) )
                  align.isvalid = true
               else
                  align = spliced_rev_extend( p, lib, geneind, curedge, read, ridx, lkmer, align ) #TODO check
                  break
               end
         elseif sg.edgetype[nodeidx] == EDGETYPE_LR || 
                sg.edgetype[nodeidx] == EDGETYPE_RR # 'LR' || 'RR'
                
               if isnull(passed_edges) # now we have to malloc
                  passed_edges = Nullable(Vector{Coordint}())
               end
               passed_extend = 0
               passed_mismat = 0
               push!(get(passed_edges), nodeidx)
               #@bp
               sgidx   -= 1
               ridx    += 1
               nodeidx -= 1
               unshift!( align.path, SGNode( geneind, nodeidx ) )
         elseif sg.edgetype[nodeidx] == EDGETYPE_SR && # 'SR'
                                ridx >= p.kmer_size
               # obligate spliced_extension
               lkmer = DNAKmer{p.kmer_size}(read.seq[(ridx-p.kmer_size+1):ridx])
               align = spliced_rev_extend( p, lib, geneind, curedge, read, ridx, lkmer, align ) ## TODO
               break
         elseif sg.edgetype[curedge] == EDGETYPE_SL ||
                sg.edgetype[curedge] == EDGETYPE_LS # 'SL' || 'LS'
               # end of alignment
               break # ?
         else #'LL' || 'RS'
               # ignore 'LL' and 'RS'
               sgidx -= 1
               ridx  += 1 # offset the lower ridx += 1
               nodeidx -= 1
               unshift!( align.path, SGNode( geneind, nodeidx ) )  
         end
         # ignore 'N'
      else 
         # mismatch
         align.mismatches += 1
         passed_mismat += 1
      end
      ridx  -= 1
      sgidx -= 1
      #print(" $(read.seq[ridx-1]),$ridx\_$(sg.seq[sgidx-1]),$sgidx ")
   end

   # if edgemat < K, spliced_extension for each in length(edges)
   # TODO go through passed_edges!!
   if !isnull(passed_edges)
      if passed_extend < p.kmer_size
         # go back.
         ext_len = Int[sg.nodelen[ i ] for i in get(passed_edges)]
         ext_len[end] = passed_extend
         rev_cumarray!(ext_len)
         #@bp
         for c in length(get(passed_edges)):-1:1  #most recent edge first
            if ext_len[c] >= p.kmer_size
               align.isvalid = true
               break
            else
               pop!( align.path ) # extension into current node failed
               align.mismatches -= passed_mismat
               cur_ridx = ridx - (sgidx - sg.nodeoffset[ get(passed_edges)[c] ]) #TODO untouched
               (cur_ridx + p.kmer_size - 1) < readlen || continue
            
               #lkmer = DNAKmer{p.kmer_size}(read.seq[(ridx-p.kmer_size):(ridx-1)])
               rkmer = DNAKmer{p.kmer_size}(read.seq[cur_ridx:(cur_ridx+p.kmer_size-1)])
               #println("$(read.seq[(ridx-p.kmer_size):(ridx-1)])\n$lkmer\n$(sg.edgeleft[get(passed_edges)[c]])")

               res_align = spliced_extend( p, lib, geneind, get(passed_edges)[c], read, cur_ridx, rkmer, align )
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

#TODO whole function needs fixing.
function spliced_rev_extend( p::AlignParam, lib::GraphLib, geneind, edgeind, read, ridx, lkmer, align )
   # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
   # returns right nodes with matching genes
   isdefined( lib.edges.left, kmer_index(lkmer) ) || return DEF_ALIGN

   best = DEF_ALIGN
   right_kmer = lib.graphs[geneind].edgeright[edgeind] # TODO
   left_nodes = intersect_sorted( lib.edges.left[ kmer_index(left_kmer) ],
                                  lib.edges.right[ kmer_index(rkmer) ], right=false)

   # do a test for trans splicing, and reset right_nodes
   for rn in left_nodes
      rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node] # + nodelen 
      res_align = ungapped_rev_extend( p, lib, rn.gene, Int(rn_offset), read, ridx, 
                                       align=deepcopy(align), nodeidx=rn.node )
      best = res_align > best ? res_align : best
   end
   if score(best) >= score(align) + p.kmer_size
      best.isvalid = true
   end
   best
end
