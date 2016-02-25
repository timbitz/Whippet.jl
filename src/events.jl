
bitstype 8 EdgeMotif

Base.convert( ::Type{EdgeMotif}, motif::UInt8 ) = box(EdgeMotif, unbox(UInt8, motif ))
Base.convert( ::Type{UInt8}, motif::EdgeMotif ) = box(UInt8, unbox(EdgeMotif, motif ))

# OBLIGATE MOTIFS (provided the cognate nodes both have reads):

const TXST_MOTIF = convert(EdgeMotif, 0b000)
const TXEN_MOTIF = convert(EdgeMotif, 0b001)

const ALTF_MOTIF = convert(EdgeMotif, 0b010)
const ALTL_MOTIF = convert(EdgeMotif, 0b011)

# REQUIRES SPANNING EDGE:

const RETI_MOTIF = convert(EdgeMotif, 0b100) 
const SKIP_MOTIF = convert(EdgeMotif, 0b101)

const ALTD_MOTIF = convert(EdgeMotif, 0b110)
const ALTA_MOTIF = convert(EdgeMotif, 0b111)

# NO MOTIF
const NONE_MOTIF = convert(EdgeMotif, 0b1000)

const MOTIF_TABLE = fill(NONE_MOTIF, 2^6 )
      # Alt TxStart  SL SL
      MOTIF_TABLE[ 0b000000 + 1 ] = TXST_MOTIF
      # Alt PolyA    RS RS
      MOTIF_TABLE[ 0b011011 + 1 ] = TXEN_MOTIF

      # Alt FirstEx  SL LS
      MOTIF_TABLE[ 0b000010 + 1 ] = ALTF_MOTIF
      # Alt LastEx   SR RS
      MOTIF_TABLE[ 0b001011 + 1 ] = ALTL_MOTIF

      # RetainedInt  LL RR
      MOTIF_TABLE[ 0b101110 + 1 ] = RETI_MOTIF

      # Skipped      RR LL
      MOTIF_TABLE[ 0b110101 + 1 ] = SKIP_MOTIF
      #              RR LR
      MOTIF_TABLE[ 0b110100 + 1 ] = SKIP_MOTIF
      #              LR LL
      MOTIF_TABLE[ 0b100101 + 1 ] = SKIP_MOTIF
      #              LR LR
      MOTIF_TABLE[ 0b100100 + 1 ] = SKIP_MOTIF

      # Alt Donor    LL LL
      MOTIF_TABLE[ 0b101101 + 1 ] = ALTD_MOTIF
      #              LL LR
      MOTIF_TABLE[ 0b101100 + 1 ] = ALTD_MOTIF

      # Alt Acceptor RR RR
      MOTIF_TABLE[ 0b110110 + 1 ] = ALTA_MOTIF
      #              LR RR
      MOTIF_TABLE[ 0b100110 + 1 ] = ALTA_MOTIF

Base.convert(::Type{EdgeMotif}, tup::Tuple{EdgeType,EdgeType}) = Base.convert(EdgeMotif, tup... )
Base.convert(::Type{EdgeMotif}, current::EdgeType, next::EdgeType) = MOTIF_TABLE[ (UInt8(current) << 3) | UInt8(next) + 1 ]

isobligate(  motif::EdgeMotif ) = motif != NONE_MOTIF && !( UInt8(motif) & 0b100 == 0b100 )
isaltsplice( motif::EdgeMotif ) = (UInt8(motif) & 0b110) == 0b110 

function process_sgquant( lib::GraphLib, graphq::GraphLibQuant )
   for g in 1:length(lib.graphs)
      process_events!( lib.graphs[g], graphq.quant[g] )
   end
end

isspanning{I <: AbstractInterval}( edge::I, node::Coordint ) = edge.first < node < edge.last ? true : false
isconnecting{I <: AbstractInterval}( edge::I, node::Coordint ) = edge.first == node || edge.last == node ? true : false

sorted_in

# this is meant for short arrays when it is faster
# than using the overhead of a set
# DEPRECATED, IntSet is very efficient
function unique_push!{T}( arr::Vector{T}, el::T )
   if !( el in arr )
      push!( arr, el )
   end
end

function hasintersect( a::IntSet, b::IntSet )
   seta,setb = length(a) > length(b) ? (b,a) : (a,b)
   for elem in seta
      if elem in setb
         return true
      end
   end
   false
end

# This function seeks to join graphs that
# should not be disjoint
function reduce_graph!( vec::Vector{IntSet}, cnt::Vector{Float64}, len::Vector{Float64}  )
   i = 1
   while i < length(vec)
      if hasintersect( vec[i], vec[i+1] )
         vec[i+1] = union( vec[i], vec[i+1] )
         cnt[i+1] += cnt[i]
         len[i+1] += len[i]
         shift!(vec)
         shift!(cnt)
         shift!(len)
         i -= 1
      end
      i += 1
   end
   vec
end

function _process_spliced( sg::SpliceGraph, sgquant::SpliceGraphQuant, node::Coordint, motif::EdgeMotif, eff_len::Int )
   inc_cnt = 0.0
   inc_len = 0.0
   inc_set = IntSet()
   exc_cnt = Vector{Float64}()
   exc_len = Vector{Float64}()
   exc_set = Vector{IntSet}()

   # push initial graph structure for inclusion/exclusion-set
   for edg in intersect( sgquant.edge, (node, node) )
      if   isconnecting( edg, node )
         inc_cnt += edg.value
         inc_len += 1
         push!( inc_set, edg.first )
         push!( inc_set, edg.last  )
      elseif isspanning( edg, node )
         push!( exc_cnt, edg.value )
         push!( exc_len, 1.0 )
         push!( exc_set, IntSet([edg.first, edg.last]) )
         reduce_graph!( exc_set, exc_cnt, exc_len )
      else
         error("Edge has to be connecting or spanning!!!!" )
      end
   end

   # if the min or max of any exclusion set is greater than the min/max
   # of the inclusion set we have a disjoint graph module and we can go
   # ahead and try to bridge nodes by iteratively adding ambiguous edges
    

   # lets finish up now.
   if exc_len == 0 # no spanning edge
      # check if we have both inclusion edges represented, or one if alt 5'/3'
      if (inc_len >= 2 && inc_cnt >= 1) || (inc_len >= 1 && inc_cnt >= 1  && isaltsplice(motif)) 
         # psi = 0.99 && likelihood_ci( psi, inc_cnt, z=1.64 )
      else
         # NA is ignored.
      end
   else # do EM
      psi = rec_spliced_em!(  )
   end
   #conf_int = likelihood_ci()
end


function output_psi{F <: AbstractFloat}( icnt::F, ecnt::F, ilen::F, elen::F, 
                                         nodes::Vector{Coordint}, node::Coordint )
   
end

function process_events( sg::SpliceGraph, sgquant::SpliceGraphQuant, eff_len::Int )
   # Heres the plan:
   # step through sets of edges, look for edge motifs, some are obligate calculations
   # others we only calculate psi if there is an alternative edge
   # Then calculate bias of the event 
   if length( sg.edges ) <= 2
      return 0
   end
   for i in 1:length(sg.edges)-1
      motif = convert(EdgeMotif, sg.edges[i], sg.edges[i+1] )
      motif == NONE_MOTIF && continue
      if isobligate( motif ) # is utr event
          
      else  # is a spliced node
         _process_spliced( sg, sgquant, sg.edges[i], motif, eff_len )
      end  
   end
end


