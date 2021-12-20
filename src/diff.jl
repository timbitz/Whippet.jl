
struct PosteriorPsi
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
function PosteriorPsi( set::Vector{PosteriorPsi}; 
                       paired::Bool=false, 
		       size::Int=1000,
                       point_est::Bool=true, 
		       pseudo_adj::Float64=0.025,
		       max_variance::Float64=0.15)
   if point_est
      psiarr = [x.psi for x in set]
      if length(set) > 1 && var(psiarr) > 0.0
         emperical = psiarr
      elseif length(set) == 1
         one = psiarr[1]
         two = one < (1.0 - pseudo_adj) ? one + pseudo_adj : one - pseudo_adj
         emperical = vcat(psiarr, two)
      else
         for i in 2:length(psiarr)
	    rand_adj = (rand() * 2pseudo_adj - pseudo_adj)
            if 0.0 > rand_adj + psiarr[i] ||
	       1.0 < rand_adj + psiarr[i]
	       rand_adj = rand_adj * -1.0
	    end
	    psiarr[i] += rand_adj
	 end
	 emperical = psiarr
      end
   else
      emperical = vcat( map( x->x.emperical, set )... )
   end
   
   if var(emperical) >= max_variance
      beta = Beta(0.1,0.1)
   else
      try
         beta = fit(Beta, emperical)
      catch e
         println(stderr, "Error: Couldn't fit beta pdf to discordant psi values:")
         println(stderr, emperical)
	 println(stderr, "Decrease `--max-variance`.")
         error(e)
      end
   end
   psi  = mean(beta)
   if !paired || point_est
      emperical = rand(beta, size)
   end
   PosteriorPsi( beta, emperical, psi )
end

# Calculates the probability, P( a-b > amt)
function probability( a::PosteriorPsi, b::PosteriorPsi; amt=0.0, iter=10 )
   success = 0
   total = min(length(a.emperical),length(b.emperical))
   for i in 1:iter
      shuffle!(a.emperical)
      shuffle!(b.emperical)
      for j in 1:total
         success += a.emperical[j] - b.emperical[j] > amt ? 1 : 0
      end
   end
   success / (total * iter)
end

function open_stream( filename::String; bufsize::Int=4_194_304 )
   fopen = open( filename, "r" )
   if isgzipped( filename )
      stream = ZlibInflateInputStream( fopen, reset_on_end=true, bufsize=bufsize )
   else
      stream = BufferedStreams.BufferedInputStream( fopen, bufsize * 4 )
   end
   stream
end

function open_streams( files::Vector{String}; bufsize::Int=4_194_304 )
   buf = Vector{BufferedStreams.BufferedInputStream}(undef, length(files))
   for i in 1:length(files)
      buf[i] = open_stream( files[i], bufsize=bufsize )
   end
   buf
end

parse_float_omit_text( str::S, header::String ) where {S <: AbstractString} = str != "NA" && str != header ? parse(Float64, str) : 0.0

function parse_psi_line( line::String; min_num=5, size=1000 )
   res  = split( line, '\t' )
   psi  = parse_float_omit_text( res[6], "Psi" )
   num  = parse_float_omit_text( res[9], "Total_Reads" )
   if psi < 0 || num < min_num
      post = PosteriorPsi()
      bool = false
   else
      post = PosteriorPsi( psi, num, size=size )
      bool = true
   end
   res,post,bool,num
end

mutable struct PosteriorEvent
   event::Vector{SubString{String}}
   complexity::String
   a::PosteriorPsi
   b::PosteriorPsi
end

parse_complexity( c::S ) where {S <: AbstractString} = split( c, COMPLEX_CHAR, keepempty=false )[1] |> x->parse(Int,x)

function process_psi_line( streams::Vector{BufferedStreams.BufferedInputStream}; 
                           min_reads=5,
			   size=1000,
			   whitelist_events = String["TS", "TE", "RI", "CE", "AA", "AD"])
   postvec = Vector{PosteriorPsi}()
   event   = split( "", "" )
   complex = 0
   entropy = 0.0
   i = 1
   mean_num = 0.0
   curnode = ""
   while i <= length(streams)
      line = readline( streams[i] )
      i += 1
      if line != ""
         par,post,isok,num = parse_psi_line( line, min_num=min_reads, size=size )
         event = par[1:5]
         event[5] in whitelist_events || (i -= 1; continue)
         !isok && continue
	 if curnode == ""
	    curnode = event[2]
	 else
	    @assert event[2] == curnode
	 end
	 mean_num += num
         push!( postvec, post )
         parcomplex = parse_complexity( par[10] )
         parentropy = parse_float_omit_text( par[11], "Entropy" )
         complex = parcomplex > complex ? parcomplex : complex
         entropy = parentropy > entropy ? parentropy : entropy
      end
   end
   mean_num /= length(postvec)
   event,complex,entropy,mean_num,postvec
end

function process_psi_files( outfile, a::Vector{BufferedStreams.BufferedInputStream}, 
                                     b::Vector{BufferedStreams.BufferedInputStream}; 
                                     min_samp::Int=1, 
				     min_reads::Int=5, 
                                     amt::Float64=0.0, 
				     size::Int=1000,
                                     point_est::Bool=true, 
				     pseudo_adj::Float64=0.01,
				     max_variance::Float64=0.15 )
   io = open( outfile, "w" )
   stream = ZlibDeflateOutputStream( io )
   output_diff_header( stream )
   while true # go through all lines until we hit eof
      a_event,a_complex,a_entropy,a_readnum,a_postvec = process_psi_line( a, min_reads=min_reads, size=size )
      b_event,b_complex,b_entropy,b_readnum,b_postvec = process_psi_line( b, min_reads=min_reads, size=size )
      if a_event == split( "", "" ) || b_event == split( "", "" )
         break # eof
      end
      complex = a_complex > b_complex ? a_complex : b_complex
      entropy = a_entropy > b_entropy ? a_entropy : b_entropy
      @assert( a_event == b_event, "Incorrect events matched!!" )
      if length(a_postvec) >= min_samp && length(b_postvec) >= min_samp
         a_post   = PosteriorPsi( a_postvec, point_est=point_est, pseudo_adj=pseudo_adj, max_variance=max_variance )
         b_post   = PosteriorPsi( b_postvec, point_est=point_est, pseudo_adj=pseudo_adj, max_variance=max_variance ) 
         fwdprob  = probability( a_post, b_post, amt=amt )
         revprob  = probability( b_post, a_post, amt=amt )
         prob     = fwdprob > revprob ? fwdprob : revprob
         psi_a    = mean(a_post.beta)
         psi_b    = mean(b_post.beta)
         deltapsi = psi_a - psi_b
         output_diff( stream, 
	              a_event, 
		      complex, 
		      entropy,
		      length(a_postvec),
		      length(b_postvec),
		      a_readnum,
		      b_readnum,
		      psi_a, 
		      psi_b, 
		      deltapsi, 
		      prob )
      end
   end
   close(stream)
   close(io)
end
