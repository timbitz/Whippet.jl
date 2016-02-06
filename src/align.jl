# requires
# using FMIndexes
# include("graph.jl")
# include("index.jl")

import Base.<,Base.>,Base.<=,Base.>=

immutable AlignParam
   mismatches::Int    # Allowable mismatches
   kmer_size::Int  # Minimum number of matches to extend past an edge
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

abstract UngappedAlignment

type SGAlignment <: UngappedAlignment
   matches::Int
   mismatches::Int
   offset::Int
   path::Vector{SGNodeTup}
   isvalid::Bool
end

const DEF_ALIGN = SGAlignment(0, 0, 0, SGNodeTup[], false)

score( align::SGAlignment ) = align.matches - align.mismatches 

>( a::SGAlignment, b::SGAlignment ) = >( score(a), score(b) )
<( a::SGAlignment, b::SGAlignment ) = <( score(a), score(b) )
>=( a::SGAlignment, b::SGAlignment ) = >=( score(a), score(b) )
<=( a::SGAlignment, b::SGAlignment ) = <=( score(a), score(b) )

function fetch_seed( p::AlignParam, read::SeqRecord, offset=p.seed_buffer )
   seed = read.seq[offset:(offset+p.seed_length-1)]
   seed
end

function try_seed( p::AlignParam, fm::FMIndex, read::SeqRecord )
   sa = 1:0
   cnt = Inf
   ctry = 0
   curpos = p.seed_buffer
   while( cnt > p.seed_tolerate && ctry <= p.seed_try && curpos <= p.seed_maxoff )
      query,offset = fetch_seed( p, read, curpos )
      sa = FMIndexes.sa_range( query, index )
      cnt = length(sa)
      ctry += 1
      curpos += p.seed_inc
   end
   sa,curpos
end

function ungapped_align( p::AlignParam, lib::GraphLib, read::SeqRecord )
   seed,readloc = try_seed( p, lib.index, read )
   locit = FMIndexes.LocationIterator( seed )
   res   = SGAlignment[]
   for s in locit
      geneind = offset_to_name( lib, s )
      align = ungapped_fwd_extend( p, lib, geneind, s - lib.offset[geneind], 
                                   read, readloc + p.seed_length - 1 ) # TODO check
      if align.isvalid
         push!(res, align)
      end
   end
   res
end


# This is the main ungapped alignment extension function in the --> direction
# Returns: SGAlignment
function ungapped_fwd_extend( p::AlignParam, lib::GraphLib, geneind, sgoffset::Int, 
                              read::SeqRecord, readoffset::Int;
                              align=SGAlignment(p.seed_length,0,sgoffset,SGNodeTup[],false),
                              nodeidx=search_sorted(lib.graphs[geneind].nodeoffset,sgoffset,lower=true) )
   ridx     = readoffset
   sgidx    = sgoffset
   sg       = lib.graphs[geneind]
   readlen  = length(read.seq)
   curnode  = nodeidx

   push!( align.path, SGNodeTup( geneind, curnode ) ) # starting node

   while( mis <= p.mismatches && ridx <= readlen )
      if read[ridx] == sg.seq[sgidx]
         # match
         align.matches += 1
      elseif (UInt8(sg.seq[sgidx]) & 0b100) == 0b100 # N,L,R,S
         if     sg.edgetype[curnode+1] == EDGETYPE_LR && 
                sg.nodelen[curnode+1]  >= p.kmer_size && # 'LR' && nodelen >= K
                readlen - ridx + 1     >= p.kmer_size
               # check edgeright[curnode+1] == read[ nextkmer ]
               # move forward K, continue 
               # This is a std exon-exon junction, if the read is an inclusion read
               # then we simply jump ahead, otherwise lets try to find a compatible right
               # node
               rkmer = SGKmer{p.kmer.size}(read[ridx:(ridx+p.kmer.size-1)])
               if lib.edges.edgeright[curnode+1] == rkmer
                  align.matches += p.kmer.size
                  ridx += 1 + p.kmer.size
                  sidx += 1 + p.kmer.size
                  curnode += 1
                  push!( align.path, SGNodeTup( geneind, curnode ) )
               else
                  align = spliced_extend( p, lib, geneind, curnode+1, read, ridx, rkmer, align )
                  break
               end
         elseif sg.edgetype[curedge] == EDGETYPE_LR || 
                sg.edgetype[curedge] == EDGETYPE_LL # 'LR' || 'LL'
               # if length(edges) > 0, push! edgematch then set to 0
               # push! edge,  (here we try to extend fwd, but keep track
               # of potential edges we pass along the way.  When the alignment
               # dies if we have not had sufficient matches we then explore
               # those edges. 
                
         elseif sg.edgetype[curedge] == EDGETYPE_LS # 'LS'
               # obligate spliced_extension
               rkmer = SGKmer{p.kmer.size}(read[ridx:(ridx+p.kmer.size-1)])
               align = spliced_extend( p, lib, geneind, curnode+1, read, ridx, kmer, align )
               break
         elseif sg.edgetype[curedge] == EDGETYPE_SR ||
                sg.edgetype[curedge] == EDGETYPE_RS # 'SR' || 'RS'
               # end of alignment
         elseif !(sg.seq[sgind] == SG_N) #'RR' || 'SL'
               # ignore 'RR' and 'SL'   
         end
         # ignore 'N'
      else 
         # mismatch
         mis += 1
      end
      ridx += 1
      sgidx += 1
   end

   # if edgemat < K, spliced_extension for each in length(edges)

   if score(align) >= p.score_min && length(align.path) == 1
      align.isvalid = true
   end
   align
end

function spliced_extend( p::AlignParam, lib::GraphLib, geneind, edgeind, read, ridx, rkmer, align )
   # Choose extending node through intersection of lib.edges.left ∩ lib.edges.right
   right_nodes = lib.graphs[geneind].edgeleft[edgeind] ∩ rkmer # returns right nodes with matching genes
   best = DEF_ALIGN
   # do a test for trans splicing, and reset right_nodes
   for rn in right_nodes
      rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node]
      res_align = ungapped_fwd_extend( p, lib, rn.gene, rn_offset, read, ridx, 
                                               align=deep_copy(align), nodeidx=rn.node )
      best = res_align > best ? res_align : best
   end
   if score(best) >= score(align)+p.kmer_size
      best.isvalid = true
   end
   best
end
