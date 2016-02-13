
using IntervalTrees

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


function increment!( quant::SpliceGraphQuant, align::Nullable{SGAlignment}; val=1.0 )
   if !isnull(align)
      if length(get(align).path) == 1
         # access node -> [ SGNode( gene, *node* ) ]
         quant.node[ get(align).path[1][2] ] += val
      else
         for n in 1:get(align).path
            
         end
      end
   end
end
