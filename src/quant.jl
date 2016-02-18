
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

function calculate_tpm!( quant::GraphLibQuant, counts::Vector{Float64}=quant.count; readlen=50, sig=0 )
   for i in 1:length(counts)
      quant.tpm[ i ] = counts[i] * readlen / quant.length[i]
   end
   const rpk_sum = sum( quant.tpm )
   for i in 1:length(quant.tpm)
      if sig > 0
         quant.tpm[i] = signif( quant.tpm[i] * SCALING_FACTOR / rpk_sum, sig )
      else
         quant.tpm[i] = ( quant.tpm[i] * SCALING_FACTOR / rpk_sum )
      end
   end
end

type Multimap
   align::Vector{SGAlignment}
   prop::Vector{Float64}
   prop_sum::Float64
end

Multimap( align::SGAlignment ) = length(align.path) >= 1 ? 
                                    Multimap( align, ones(length(align.path)) / length(align.path) ) :
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

function unsafe_copy!{T <: Number}( dest::Vector{T}, src::Vector{T} )
   for i in 1:length(src)
      dest[i] = src[i]
   end
end 

# This function recursively performs expectation maximization
# and sets the quant.tpm array as the proposed expression set
# at each iteration.   It also requires that calculate_tpm! be 
# run initially once prior to rec_gene_em!() call.
@fastmath function rec_gene_em!( quant::GraphLibQuant,
                                 ambig::Vector{Multimap}, 
                                 temp::Vector{Float64}=ones(length(quant.count)),
                                 uniqsum=sum(quant.count),
                                 ambigsum=length(ambig),
                                 it=1, max=1000, sig=0)
   
   unsafe_copy!( temp, quant.count )
   for mm in ambig
      if it > 1 # Maximization
         mm.prop_sum = 0.0
         for ai in 1:length(mm.align)
            init_gene = mm.align[ai].path[1].gene
            init_tpm  = quant.tpm[ init_gene ]
            mm.prop[ai] = init_tpm
            mm.prop_sum += init_tpm
         end
      end
      
      for ai in 1:length(mm.align)
         init_gene = mm.align[ai].path[1].gene
         temp[ init_gene ] += (mm.prop[ai] / mm.prop_sum)
      end
   end

   calculate_tpm!( quant, temp, sig=sig ) # Expectation
     
   #iterate
   if temp != quant.tpm && it < max
      it = rec_gene_em( quant, ambig, temp, uniqsum=uniqsum, ambigsum=ambigsum, 
                        it=it+1, max=max, sig=sig)
   end
   it
end


