


immutable AlignParam
   mismatches::Int # Allowable mismatches
   seednum::Int    # Starting number of seeds
   seedtol::Int    # Seed matching tolerance 
   seedlen::Int    # Seed Size
   seedreg::Int    # Seed gathering constraint
   stranded::Bool  # Is input data + strand only?
   paired::Bool    # Paired end data?
end


