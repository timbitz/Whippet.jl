# requires
# using FMIndexes
# include("graph.jl")
# include("index.jl")

import Base.<,Base.>,Base.<=,Base.>=

import Bio.Seq.SeqRecord

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

AlignParam() = AlignParam( 2, 9, 3, 7, 18, 35, 5, 10, 1, 10, 25, true, false, false, true )

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


@debug function try_seed( p::AlignParam, index::FMIndex, read::SeqRecord )
   sa = 2:1
   cnt = 1000
   ctry = 0
   curpos = p.seed_buffer
   while( cnt > p.seed_tolerate && ctry <= p.seed_try && curpos <= p.seed_maxoff )
      sa = FMIndexes.sa_range( read.seq[curpos:(curpos+p.seed_length-1)], index )
      @bp
      cnt = length(sa)
      ctry += 1
      curpos += p.seed_inc
   end
   sa,curpos
end

@debug function ungapped_align( p::AlignParam, lib::GraphLib, read::SeqRecord )
   seed,readloc = try_seed( p, lib.index, read )
   locit = FMIndexes.LocationIterator( seed, lib.index )
   res   = Nullable{Vector{SGAlignment}}()
   for s in locit
      geneind = search_sorted( lib.offset, convert(Coordint, s), lower=true ) 
      @bp
      align = ungapped_fwd_extend( p, lib, geneind, s - lib.offset[geneind], 
                                   read, readloc + p.seed_length ) # TODO check
      if align.isvalid
         if isnull( res )
            res = Nullable(Vector{SGAlignment})
         end
         push!(res, align)
      end
   end
   res
end


# This is the main ungapped alignment extension function in the --> direction
# Returns: SGAlignment
@debug function ungapped_fwd_extend( p::AlignParam, lib::GraphLib, geneind, sgidx::Int, 
                              read::SeqRecord, ridx::Int;
                              align=SGAlignment(p.seed_length,0,sgidx,SGNodeTup[],false),
                              nodeidx=search_sorted(lib.graphs[geneind].nodeoffset,Coordint(sgidx),lower=true) )
   sg       = lib.graphs[geneind]
   readlen  = length(read.seq)
   curedge  = nodeidx

   passed_edges = Nullable{Vector{Coordint}}() # don't allocate arrays unless
   edge_matches = Nullable{Vector{Coordint}}() # we are going to need them
   extend_match = 0

   push!( align.path, SGNodeTup( geneind, nodeidx ) ) # starting node

   while( align.mismatches <= p.mismatches && ridx <= readlen )
      @bp
      if read.seq[ridx] == sg.seq[sgidx]
         # match
         align.matches += 1
         extend_match  += 1
      elseif (UInt8(sg.seq[sgidx]) & 0b100) == 0b100 # N,L,R,S
         curedge = nodeidx+1
         if     sg.edgetype[curedge] == EDGETYPE_LR && 
                sg.nodelen[curedge]  >= p.kmer_size && # 'LR' && nodelen >= K
                readlen - ridx + 1   >= p.kmer_size
               # check edgeright[curnode+1] == read[ nextkmer ]
               # move forward K, continue 
               # This is a std exon-exon junction, if the read is an inclusion read
               # then we simply jump ahead, otherwise lets try to find a compatible right
               # node
               rkmer = SGKmer{p.kmer.size}(read.seq[ridx:(ridx+p.kmer.size-1)])
               if lib.edges.edgeright[curedge] == rkmer
                  align.matches += p.kmer.size
                  ridx += p.kmer.size
                  sidx += 1 + p.kmer.size
                  nodeidx += 1
                  push!( align.path, SGNodeTup( geneind, nodeidx ) )
               else
                  align = spliced_extend( p, lib, geneind, curedge, read, ridx, rkmer, align )
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
                  edge_matches = Nullable(Vector{Coordint}())
                  extend_match = 0
               else
                  push!(get(edge_matches), extend_match)
               end
               push!(get(passed_edges), curedge)
               extend_match = 0
               sidx += 1
               ridx -= 1
         elseif sg.edgetype[curedge] == EDGETYPE_LS # 'LS'
               # obligate spliced_extension
               rkmer = SGKmer{p.kmer.size}(read.seq[ridx:(ridx+p.kmer.size-1)])
               align = spliced_extend( p, lib, geneind, curedge, read, ridx, kmer, align )
               break
         elseif sg.edgetype[curedge] == EDGETYPE_SR ||
                sg.edgetype[curedge] == EDGETYPE_RS # 'SR' || 'RS'
               # end of alignment
               break # ?
         elseif !(sg.seq[sgind] == SG_N) #'RR' || 'SL'
               # ignore 'RR' and 'SL'
               sidx += 1
               ridx -= 1 # offset the lower ridx += 1  
         end
         # ignore 'N'
      else 
         # mismatch
         align.mismatches += 1
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
   if score(best) >= score(align) + p.kmer_size
      best.isvalid = true
   end
   best
end
