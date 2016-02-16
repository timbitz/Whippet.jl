
using IntervalTrees

const SCALING_FACTOR = 1_000_000

# This is where we count reads for nodes/edges/circular-edges
immutable SpliceGraphQuant
   node::Vector{Float64}
   edge::IntervalMap{Exonmax,Float64}
   circ::Dict{Tuple{Exonmax,Exonmax},Float64}
end

# Default constructer
SpliceGraphQuant() = SpliceGraphQuant( Vector{Float64}(),
                                       IntervalMap{Exonmax,Float64}(),
                                       Dict{Tuple{Exonmax,Exonmax},Float64}() )

SpliceGraphQuant( sg::SpliceGraph ) = SpliceGraphQuant( zeros( length(sg.nodelen) ),
                                                        IntervalMap{Exonmax,Float64}(),
                                                        Dict{Tuple{Exonmax,Exonmax},Float64}() )

function increment!( quant::SpliceGraphQuant, align::SGAlignment; val=1.0 )
   align.isvalid == true || return
   if length(align.path) == 1
      # access node -> [ SGNode( gene, *node* ) ]
      quant.node[ align.path[1].node ] += val
   else
      # TODO: potentially handle trans-splicing here via align.istrans variable?

      # Otherwise, lets step through pairs of nodes and add val to those edges
      for n in 1:(length(align.path)-1)
         lnode = align.path[n].node
         rnode = align.path[n+1].node
         if lnode < rnode
            interv = Interval{Exonmax}( lnode, rnode )
            quant.edge[ interv ] = get( quant.edge, interv, IntervalValue(0,0,0.0) ).value + val
         elseif lnode >= rnode
            quant.circ[ (lnode, rnode) ] = get( quant.circ, (lnode,rnode), 0.0) + val
         end
      end
   end
end


function transcripts_per_mil!( tpm::Vector{Float64}, counts::Vector{Float64}, lengths::Vector{Float64}; 
                               readlen=50 )
   for i in 1:length(counts)
      tpm[ i ] = counts[i] * readlen / lengths[i]
   end
   const rpk_sum = sum( tpm )
   for i in 1:length(tpm)
      tpm[i] = ( tpm[i] * SCALING_FACTOR / rpk_sum )
   end
end


