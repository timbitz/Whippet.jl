
immutable PosteriorPsi
   beta::Distributions.Beta
   emperical::Vector{Float64}
   psi::Float64
end

# create one PosteriorPsi object given psi, N, and priors α and β
function PosteriorPsi( psi::Float64, N::Float64, α=1.0, β=1.0; size=1000 )
   # Posterior Beta is k + α, n-k + β
   beta = Beta( psi*N + α, (1-psi)*N + β )
   emperical = rand(beta, size)
   psi = mean(beta)
   PosteriorPsi( beta, emperical, psi )
end

# combine multiple beta posteriors for replicates into one
function PosteriorPsi( set::Vector{PosteriorPsi} )
   emperical = vcat( map( x->x.emperical, set )... )
   beta = fit(Beta, emperical)
   psi = mean(beta)
   PosteriorPsi( beta, emperical, psi )
end

# Calculates the probability, P( a-b > amt)
function probability( a::PosteriorPsi, b::PosteriorPsi; amt=0.0 )
   good = 0
   total = min(length(a.emperical),length(b.emperical))
   for i in 1:total
      good += a.emperical[i] - b.emperical[i] > amt ? 1 : 0
   end
   @fastmath good / total
end

function open_stream( filename )
   fopen = open( filename, "r" )
   if isgzipped( filename )
      stream = ZlibInflateInputStream( fopen, reset_on_end=true )
   else
      stream = BufferedStreams.BufferedInputStream( fopen )
   end
   stream
end

function open_streams( files::Vector{ASCIIString} )
   buf = Vector{BufferedStreams.BufferedInputStream}(length(files))
   for i in 1:length(files)
      buf[i] = open_stream( files[i] )
   end
end

function parse_psi_line( line::ASCIIString )
   res = split( line, '\t' )
   
end

function process_psi_line( streams::Vector{BufferedStreams.BufferedInputStream} )
   
   for bs in streams
      line = readline( bs )
      if line != ""
         
      end
   end
end

function process_psi_files( a::Vector{BufferedStreams.BufferedInputStream}, 
                            b::Vector{BufferedStreams.BufferedInputStream} )
   i = 0
   
end
