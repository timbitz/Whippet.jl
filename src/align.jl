
# requires
# using FMIndexes
# include("graph.jl")
# include("index.jl")

immutable AlignParam
   mismatches::Int    # Allowable mismatches
   seed_try::Int      # Starting number of seeds
   seed_tolerate::Int # Allow at most _ hits for a valid seed
   seed_length::Int   # Seed Size
   seed_maxoff::Int   # Seed gathering constraint
   seed_buffer::Int   # Ignore first _ bases
   seed_inc::Int      # Incrementation for subsequent seed searches
   seed_rng::Int      # Seed for random number generator
   score_range::Int   # Scoring range to be considered repetitive
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

function try_seed( fm::FMIndex, p::AlignParam, read::SeqRecord )
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

function ungapped_align( lib::GraphLib, param::AlignParam, read::SeqRecord )
   seed,readloc = try_seed( lib.index, param, read )
   gene = sorted_getname  
end

function ungapped_fwd_extend()
end
function ungapped_rev_extend()
end
 
