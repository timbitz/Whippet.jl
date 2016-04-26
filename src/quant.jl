
using IntervalTrees

const SCALING_FACTOR = 1_000_000

# This is where we count reads for nodes/edges/circular-edges/effective_lengths
# bias is an adjusting multiplier for nodecounts to bring them to the same level
# as junction counts, which should always be at a lower level, e.g. bias < 1
type SpliceGraphQuant
   node::Vector{Float64}
   edge::IntervalMap{Exonmax,Float64}
   circ::Dict{Tuple{Exonmax,Exonmax},Float64}
   leng::Vector{Float64}
   bias::Float64
end

# Default constructer
SpliceGraphQuant() = SpliceGraphQuant( Vector{Float64}(),
                                       IntervalMap{Exonmax,Float64}(),
                                       Dict{Tuple{Exonmax,Exonmax},Float64}(),
                                       Vector{Float64}(), 1.0 )

SpliceGraphQuant( sg::SpliceGraph ) = SpliceGraphQuant( zeros( length(sg.nodelen) ),
                                                        IntervalMap{Exonmax,Float64}(),
                                                        Dict{Tuple{Exonmax,Exonmax},Float64}(),
                                                        zeros( length(sg.nodelen) ), 1.0 )


# Here we store whole graphome quantification
immutable GraphLibQuant
   tpm::Vector{Float64}
   count::Vector{Float64}
   length::Vector{Float64}
   quant::Vector{SpliceGraphQuant}
end

function GraphLibQuant( lib::GraphLib, ref::Refset )
   tpm    = zeros( length(lib.graphs) )
   count  = zeros( length(lib.graphs) )
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
      const denom = max( 1.0, (quant.length[i] - readlen) )
      @fastmath quant.tpm[ i ] = counts[i] / denom
   end
   const rpk_sum = sum( quant.tpm )
   for i in 1:length(quant.tpm)
      if sig > 0
         @fastmath quant.tpm[i] = signif( quant.tpm[i] * SCALING_FACTOR / rpk_sum, sig )
      else
         @fastmath quant.tpm[i] = ( quant.tpm[i] * SCALING_FACTOR / rpk_sum )
      end
   end
end

type Multimap
   align::Vector{SGAlignment}
   prop::Vector{Float64}
   prop_sum::Float64
end

Multimap( aligns::Vector{SGAlignment} ) = length(aligns) >= 1 ? 
                                    Multimap( aligns, ones(length(aligns)) / length(aligns), 1.0 ) :
                                    Multimap( aligns, Float64[], 0.0 )


function assign_ambig!( graphq::GraphLibQuant, ambig::Vector{Multimap}; ispaired=false )
   for mm in ambig
      i = 1
      while i <= length(mm.prop)
         if ispaired
            (i == length(mm.prop)) && break
            count!( graphq, mm.align[i], mm.align[i+1], val=mm.prop[i] )
            i += 1
         else
            count!( graphq, mm.align[i], val=mm.prop[i] )
         end
         i += 1
      end
   end
end

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
         # trans-spicing off->
         align.path[n].gene != init_gene && continue
         align.path[n+1].gene != init_gene && continue
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

function calculate_bias!( sgquant::SpliceGraphQuant )
   edgecnt = 0.0
   for edgev in sgquant.edge
      edgecnt += edgev.value
   end
   @fastmath edgelevel = edgecnt / length(sgquant.edge)
   nodecnt = 0.0
   for nodev in sgquant.node
      nodecnt += nodev
   end
   @fastmath nodelevel = nodecnt / sum( sgquant.leng )
   # never down-weight junction reads, only exon-body reads
   bias = @fastmath min( edgelevel / nodelevel, 1.0 )
   sgquant.bias = bias
   bias
end

function global_bias( graphq::GraphLibQuant )
   bias_ave = 0.0
   bias_var = 0.0
   n = 1
   for sgq in graphq.quant
      curbias = calculate_bias!( sgq )
      if curbias > 10
         println( sgq )
      end
      old = bias_ave
      @fastmath bias_ave += (curbias - bias_ave) / n
      if n > 1
        @fastmath bias_var += (curbias - old)*(curbias - bias_ave)
      end
      n += 1
   end
   bias_ave,bias_var
end

# Every exon-exon junction has an effective mappable space of readlength - K*2.
# Since we only count nodes that don't map over edges, we can give
# give each junction an effective_length of one, and now the space
# inside of a node for which a read could map without overlapping a junction
# is put into units of junction derived effective_length
# kadj is minimum of (readlength - minalignlen) and k-1
@inline function eff_length( node, sg::SpliceGraph, eff_len::Int, kadj::Int )
#   len = sg.nodelen[node] + (istxstart( sg.edgetype[node] ) ? 0 : kadj) +
#                            (istxstop( sg.edgetype[node+1] ) ? 0 : kadj)
   first = sg.nodeoffset[node]
   last  = sg.nodeoffset[node] + sg.nodelen[node] - 1
   len   = sg.nodelen[node] + (istxstart( sg.edgetype[node] ) ? 0 : kadj) +
                              (istxstop( sg.edgetype[node+1] ) ? 0 : kadj)
   @fastmath map = sum(sg.map[first:last]) / length(first:last)
   @fastmath (len * map) / eff_len
end

function eff_lengths!( sg::SpliceGraph, sgquant::SpliceGraphQuant, eff_len::Int, kadj::Int )
   for i in 1:length( sg.nodelen )
      sgquant.leng[i] = eff_length( i, sg, eff_len, kadj )
   end
end

function effective_lengths!( lib::GraphLib, graphq::GraphLibQuant, eff_len::Int, kadj::Int )
   for i in 1:length( lib.graphs )
      eff_lengths!( lib.graphs[i], graphq.quant[i], eff_len, kadj )
   end
end



function Base.unsafe_copy!{T <: Number}( dest::Vector{T}, src::Vector{T} )
   for i in 1:length(src)
      dest[i] = src[i]
   end
end 

# This function recursively performs expectation maximization
# and sets the quant.tpm array as the proposed expression set
# at each iteration.   It also requires that calculate_tpm! be 
# run initially once prior to rec_gene_em!() call.
function rec_gene_em!( quant::GraphLibQuant, ambig::Vector{Multimap}; 
                                 count_temp::Vector{Float64}=ones(length(quant.count)),
                                 tpm_temp::Vector{Float64}=ones(length(quant.count)),
                                 uniqsum=sum(quant.count),
                                 ambigsum=length(ambig),
                                 it=1, max=1000, sig=0, readlen=50 )
   
   unsafe_copy!( count_temp, quant.count )
   unsafe_copy!( tpm_temp, quant.tpm )

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
         @fastmath prop = mm.prop[ai] / mm.prop_sum
         mm.prop[ai] = prop
         @fastmath count_temp[ init_gene ] += prop
      end
   end

   calculate_tpm!( quant, count_temp, sig=sig, readlen=readlen ) # Expectation
 
   #iterate
   if tpm_temp != quant.tpm && it < max
      it = rec_gene_em!( quant, ambig, count_temp=count_temp, tpm_temp=tpm_temp, uniqsum=uniqsum, ambigsum=ambigsum, 
                        it=it+1, max=max, sig=sig, readlen=readlen)
   end
   it
end

function output_tpm( file, lib::GraphLib, gquant::GraphLibQuant )
   io = open( file, "w" )
   stream = ZlibDeflateOutputStream( io )
   for i in 1:length(lib.names)
      write( stream, lib.names[i] )
      write( stream, '\t' )
      write( stream, string(gquant.tpm[i]) )
      write( stream, '\n' )
   end
   close( stream )
   close( io )
end
