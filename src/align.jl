
# requires
# using FMIndexes
# include("graph.jl")
# include("index.jl")

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
end

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
   locit = LocationIterator( seed )
   for s in locit
      geneind = offset_to_name( lib, s )
      path = ungapped_fwd_extend( lib.graphs[geneind], s - lib.offset[geneind], 
                                  read, readloc + p.seed_length - 1 ) # TODO check
   end
end

function ungapped_fwd_extend( p::AlignParam, sgarray, sgind, sgoffset::Int, 
                              read::SeqRecord, readoffset::Int )
   align    = SGAlignment( p.seed_length, 0, sgoffset, SGNodeTup[] ) 
   ridx     = readoffset
   sgidx    = sgoffset
   sg       = sgarray[sgind]
   curnode  = search_sorted( sg.nodeoffset, sgoffset, lower=true )
   curedge  = curnode+1  
 
   while( mis < p.mismatches ) # add < length(sg.seq)
      if read[ridx] == sg.seq[sgind]
         # match
         align.matches += 1
      elseif (UInt8(sg.seq[sgind]) & 0b100) == 0b100 # N,L,R,S
         if     sg.edgetype[curedge] == EDGETYPE_LR && 
                sg.nodelen[curedge] >= p.kmer_size  # 'LR' && nodelen >= K
               # check edgeright[curnode+1] == 
               # move forward K, continue    
         elseif sg.edgetype[curedge] == EDGETYPE_LR || 
                sg.edgetype[curedge] == EDGETYPE_LL # 'LR' || 'LL'
               # if length(edges) > 0, push! edgematch then set to 0
               # push! edge,       
         elseif sg.edgetype[curedge] == EDGETYPE_LS # 'LS'
               # obligate spliced_extension
         elseif sg.edgetype[curedge] == EDGETYPE_SR ||
                sg.edgetype[curedge] == EDGETYPE_RS # 'SR' || 'RS'
               # end of alignment
         elseif !(sg.seq[sgind] == SG_N) #'RR' || 'SL'
            
         end
         # ignore 'N', 'RR' and 'SL' 
      else 
         # mismatch
         mis += 1
      end
      ridx += 1
      sgidx += 1
   end

   # if edgemat < K, spliced_extension for each in length(edges)

   if (align.mat - align.min) >= p.score_min
      unshift!( align.path, SGNodeTup( sgind, curnode ) )
   else
      # return no valid alignment
   end
end

function ungapped_edge_extend( p::AlignParam, sgarray, sgnode )
   
end
