
# requires
# using FMIndexes
# include("graph.jl")
# include("index.jl")

immutable AlignParam
   mismatches::Int    # Allowable mismatches
   seed_try::Int      # Starting number of seeds
   seed_length::Int   # Seed Size
   seed_region::Int   # Seed gathering constraint
   seed_offset::Int   # Ignore first _ bases
   is_stranded::Bool  # Is input data + strand only?
   is_paired::Bool    # Paired end data?
   is_trans_ok::Bool  # Do we allow edges from one gene to another
   is_circ_ok::Bool   # Do we allow back edges
end

function check_seed_hit( index::FMIndex, param )
   sa = FMIndexes.sa_range( query, index )
   cnt = length(sa)
   FMIndexes.LocationIterator( sa, index )
end


function ungapped_align( lib::GraphLib, param::AlignParam, reads )

end

function ungapped_fwd_extend()
end
function ungapped_rev_extend()
end
 
