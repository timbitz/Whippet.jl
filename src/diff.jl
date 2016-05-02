

immutable PosteriorPsi
   beta::Distributions.Beta
   sample::Vector{Float64}
   psi::Float64
end

function PosteriorPsi( psi::Float64, N::Float64, α=1.0, β=1.0; size=1000 )
   # Posterior Beta is k + α, n-k + β
   beta = Beta( psi*N + α, (1-psi)*N + β )
   sample = rand(beta, size)
   psi = mean(beta)
   PosteriorPsi( beta, sample, psi )
end

# combine multiple beta posteriors for replicates into one
function PosteriorPsi( set::Vector{PosteriorPsi} )
   sample = vcat( map( x->x.sample, set )... )
   beta = fit(Beta, sample)
   psi = mean(beta)
   PosteriorPsi( beta, sample, psi )
end

# Calculates the probability, P(|a-b| > amt)
function probability( a::PosteriorPsi, b::PosteriorPsi; amt=0.0 )
   good = 0
   total = min(length(a.sample),length(b.sample))
   for i in 1:total
      good += abs(a.sample[i] - b.sample[i]) > amt ? 1 : 0
   end
   @fastmath good / total
end
