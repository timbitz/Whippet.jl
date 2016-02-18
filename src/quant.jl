
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


# Here we store whole graphome quantification
immutable GraphLibQuant
   tpm::Vector{Float64}
   count::Vector{Float64}
   length::Vector{Float64}
   quant::Vector{SpliceGraphQuant}
end

function GraphLibQuant( lib::GraphLib, ref::Refset )
   tpm    = zeros( length(lib.graphs) )
   count  =  ones( length(lib.graphs) )
   len    =  ones( length(lib.graphs) )
   quant  = Vector{SpliceGraphQuant}( length(lib.graphs) )
   for i in 1:length( lib.graphs )
      name = lib.names[i]
      if haskey( ref.geneset, name )
         len[i] = ref.geneset[name].length
      end
      quant[i] = SpliceGraphQuant( lib.graphs[i] )
   end
   GraphLibQuant( tpm, count, len, quant )
end

function calculate_tpm!( quant::GraphLibQuant; readlen=50 )
   for i in 1:length(quant.count)
      quant.tpm[ i ] = quant.count[i] * readlen / quant.length[i]
   end
   const rpk_sum = sum( quant.tpm )
   for i in 1:length(quant.tpm)
      quant.tpm[i] = ( quant.tpm[i] * SCALING_FACTOR / rpk_sum )
   end
end

immutable Multimap
   align::Vector{SGAlignment}
   prop::Vector{Float64}
end

Multimap( align::SGAlignment ) = length(align.path) >= 1 ? 
                                    Multimap( align, ones(length(align.path)) / length(align.path) )
                                    Multimap( align, Float64[] )


function count!( graphq::GraphLibQuant, align::SGAlignment; val=1.0 )
   align.isvalid == true || return
   init_gene = align.path[1].gene
   sgquant   = graphq.quant[ init_gene ]
   
   graphq.count[ init_gene ] += 1
   
   if length(align.path) == 1
      # access node -> [ SGNode( gene, *node* ) ]
      sgquant.node[ align.path[1].node ] += val
   else
      # TODO: potentially handle trans-splicing here via align.istrans variable?
      
      # Otherwise, lets step through pairs of nodes and add val to those edges
      for n in 1:(length(align.path)-1)
         lnode = align.path[n].node
         rnode = align.path[n+1].node
         if lnode < rnode
            interv = Interval{Exonmax}( lnode, rnode )
            sgquant.edge[ interv ] = get( sgquant.edge, interv, IntervalValue(0,0,0.0) ).value + val
         elseif lnode >= rnode
            sgquant.circ[ (lnode, rnode) ] = get( sgquant.circ, (lnode,rnode), 0.0) + val
         end
      end
   end
end
