divsum{T <: AbstractFloat}( arr::Vector{T} ) = arr / convert(T, sum(arr))

function rec_em{F <: AbstractFloat, T <: Tuple}( unique_cnts::Vector{F}, 
                                                 lengths::Vector{Int}, 
                                                 prop::Vector{F}, # proposed model 
                                                 ambig::Vector{T}; #) req args 
                                                 uniqsum=sum(unique_cnts), 
                                                 ambigsum=length(ambig), 
                                                 it=1, max=1000, sig=0)
    norm = deepcopy(unique_cnts)
    for am_set in ambig
         am_tmp = Tuple{Int,Float64}[]
         tmp_sum = 0.0

         if it == 1  # start with uniform distribution
             for i in am_set
                push!(am_tmp, (i, 1/length(am_set)))
             end
             tmp_sum = Float64(length(am_set))
         else # otherwise use previous model
             for i in am_set
                @inbounds push!(am_tmp, (i,prop[i]))
                @inbounds tmp_sum += prop[i]
             end
         end

         for (i,j) in am_tmp
             @inbounds norm[i] += (j / tmp_sum)
         end
    end
    println("Inc: $(norm[1]), Exc: $(norm[2])")

    norm /= uniqsum + ambigsum
    norm_sum = 0.0
    for i in 1:length(norm)
       @inbounds norm[i] /= lengths[i]
       @inbounds norm_sum += norm[i]
    end

    if sig > 0 #adjust to 'sign'-decimal places for convergence test
       for i in 1:length(norm)
          @inbounds norm[i] = signif( norm[i]/norm_sum, sig )
       end
    else 
       norm /= norm_sum   
    end

    if norm != prop && it < max
         # iterate
         return( rec_em( unique_cnts, lengths, norm, ambig, 
                         uniqsum=uniqsum, ambigsum=ambigsum, 
                         it=it+1, max=max, sig=sig ) )
     end

    norm, it
end


function em_test()
   # genes and their non-ambiguous counts
   gene_counts = Float64[ 1,2,3,4 ]
   gene_lengths = [2,4,8,16]
   # each ambiguous read is a tuple of possible assignments (gene indexes)
   gene_ambig = [(1,2),(1,4),(2,3),(3,4)] 
   @time rec_em( gene_counts, gene_lengths, divsum(gene_counts), gene_ambig )
   # Convergence in 30 iterations
   #  0.000116 seconds (1.86 k allocations: 64.359 KB)
   #([0.4421102706288989,0.2705410604398086,0.18391422825839362,0.10343444067289885],30)


   # For Psi it is the same thing
   # Lets try 3 exons of 'length' 20 [__C1__]--------[__A__]--------[__C2__]
   psi_counts = Float64[ 5,5 ]  #inclusion reads (EEJ+A), exclusion reads (only EEJ)
   psi_ambig = [(1,2) for i in 1:10] # add 40 flanking ambiguous reads
   psi_lengths = [ 60, 40 ]  
   @time rec_em( psi_counts, psi_lengths, divsum(psi_counts), psi_ambig )

   # Calculates 63% PSI in 123 iterations
   #julia> @time rec_em( psi_counts, psi_lengths, divsum(psi_counts), psi_ambig, max=1000 )
   #  0.001680 seconds (41.95 k allocations: 1.245 MB)
   #([0.6368956325934994,0.3631043674065007],123)

   # But we only need the first 3-4 positions, so we can cap the iterations at 10
   # and get a huge bump in performance.
   #julia> @time rec_em( psi_counts, psi_lengths, divsum(psi_counts), psi_ambig, max=10 )
   #  0.000172 seconds (3.41 k allocations: 103.656 KB)
   #([0.6342180896585914,0.36578191034140856],10)

   #  [__C1__]----[__A1__]----[__A2__]----[__C2__]

   #I  -    -    -    - 
   #  [C1] [A1] [A2] [C2] = 0.25
   #J     -    -    -    -
   #   -    -    -
   #  [C1] [A2] [C2]      = 0.25
   #      -    -    -
   #   -    -    -
   #  [C1] [A1] [C2]      = 0.25
   #      -    -    -
   #---------------------------------
   #   -    -
   #  [C1] [C2]           = 0
   #     -    -

   psi_counts = Float64[ 6, 1 ]
   psi_lengths = [ 7, 5 ]
   psi_ambig = [(1,2) for i in 1:12]
   @time rec_em( psi_counts, psi_lengths, divsum(psi_counts), psi_ambig, sig=4 )

   # Test imbalanced exon body to exon-exon junction reads
   psi_counts = Float64[ 5+5, 5 ]
   psi_lengths = [ 6, 4 ]
   psi_ambig = [(1,2) for i in 1:10]
   @time rec_em( psi_counts, psi_lengths, divsum(psi_counts), psi_ambig, sig=4 )
end

