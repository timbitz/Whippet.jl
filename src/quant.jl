
using IntervalTrees

# This is where we count reads for nodes/edges/circular-edges
immutable SpliceGraphQuant
   nodecnt::Vector{Float64}
   edgecnt::IntervalMap{Exonmax,Float64}
   circcnt::Dict{Tuple{Exonmax,Exonmax},Float64}
end

# Default constructer
SpliceGraphQuant() = SpliceGraphQuant( Vector{Float64}(), Vector{UInt32}(),
                                       IntervalMap{Exonmax,Float64}(),
                                       Dict{Tuple{Exonmax,Exonmax},Float64}() )




function add!( quant::SpliceGraphQuant, align::Nullable{SGAlignment} )
   if !isnull(align)
      for n in 1:get(align).path

      end
   end
end
