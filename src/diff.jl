
immutable PosteriorPsi
   beta::Distributions.Beta
   emperical::Vector{Float64}
   psi::Float64
end

PosteriorPsi() = PosteriorPsi( Beta(1,1), Float64[], 0.0 )

# create one PosteriorPsi object given psi, N, and priors α and β
# Posterior                      Likelihood         Prior
# Beta( Ψ*N + α, (1-Ψ)*N + β ) ∝ Binomial( Ψ, N ) * Beta( α, β )
function PosteriorPsi( Ψ::Float64, N::Float64, α=1.0, β=1.0; size=1000 )
   # Posterior Beta is k + α, n-k + β
   beta = Beta( Ψ*N + α, (1-Ψ)*N + β )
   emperical = rand(beta, size)
   psi = mean(beta)
   PosteriorPsi( beta, emperical, Ψ )
end

# combine multiple beta posteriors for replicates into one
# if unpaired, we don't need to maintain sample1_rep1 to sample2_rep1
# so we can re-sample from the new joint beta
function PosteriorPsi( set::Vector{PosteriorPsi}; paired=false, size=1000 )
   emperical = vcat( map( x->x.emperical, set )... )
   beta = fit(Beta, emperical)
   psi  = mean(beta)
   if !paired
      emperical = rand(beta, size)
   end
   PosteriorPsi( beta, emperical, psi )
end

# Calculates the probability, P( a-b > amt)
function probability( a::PosteriorPsi, b::PosteriorPsi; amt=0.0 )
   success = 0
   total = min(length(a.emperical),length(b.emperical))
   for i in 1:total
      success += a.emperical[i] - b.emperical[i] > amt ? 1 : 0
   end
   @fastmath success / total
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

function parse_psi_line( line::ASCIIString; min_num=5 )
   res  = split( line, '\t' )
   psi  = res[6] != "NA" ? parse(Float64, res[6]) : 0.0
   num  = res[9] != "NA" ? parse(Float64, res[9]) : 0.0
   if psi < 0 || num <= min_num
      post = PosteriorPsi()
      bool = false
   else
      post = PosteriorPsi( psi, num )
      bool = true
   end
   res[1:5],post,bool
end

type PosteriorEvent
   event::Vector{SubString{ASCIIString}}
   complexity::ASCIIString
   a::PosteriorPsi
   b::PosteriorPsi
end

parse_complexity( c::ASCIIString ) = split( c, 'C', keep=false )[1] |> x->parse(Int,x)

function process_psi_line( streams::Vector{BufferedStreams.BufferedInputStream} )
   postvec = Vector{PosteriorPsi}()
   event   = split( "", "" )
   complex = 0
   for bs in streams
      line = readline( bs )
      if line != ""
         par,post,isok = parse_psi_line( line )
         !isok && continue
         push!( postvec, post )
         event = par[1:4]
         parcomplex = parse_complexity( par[5] )
         complex = parcomplex > complex ? parcomplex : complex
      end
   end
   event,complex,postvec
end

function process_psi_files( a::Vector{BufferedStreams.BufferedInputStream}, 
                            b::Vector{BufferedStreams.BufferedInputStream}; min_samp=1 )
   i = 0
   
   a_event,a_complex,a_post = process_psi_line( a )
   b_event,b_complex,b_post = process_psi_line( b )
   if length(a_post) > min_samp && length(b_post) > min_samp
      a_post = PosteriorPsi( a_post )
      b_post = PosteriorPsi( b_post )
   end
end
