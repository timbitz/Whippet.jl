
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

# this is meant for short arrays when it is faster
# than using the overhead of a set
function unique_push!{T}( arr::Vector{T}, el::T )
   if !( el in arr )
      push!( arr, el )
   end
end

function _process_spliced( sg::SpliceGraph, sgquant::SpliceGraphQuant, node::Coordint, motif::EdgeMotif, eff_len::Int )
   inc_cnt = 0.0
   exc_cnt = 0.0
   inc_len = 0.0
   exc_len = 0.0
   nodeset = Vector{Coordint}()

   for edg in intersect( sgquant.edge, (node, node) )
      if isconnecting( edg, node )
         inc_cnt += edg.value
         inc_len += 1
         node == edg.first || unique_push!( nodeset, edg.first )
         node == edg.last  || unique_push!( nodeset, edg.last  )
      elseif isspanning( edg, node )
         exc_cnt += edg.value
         exc_len += 1
         unique_push!( nodeset, edg.first )
         unique_push!( nodeset, edg.last  )
      else
         error("Edge has to be connecting or spanning!!!!" )
      end
   end
   
   if exc_len == 0 # no spanning edge
      # check if we have both inclusion edges represented, or one if alt 5'/3'
      if (inc_len >= 2 && inc_cnt >= 1) || (inc_len >= 1 && inc_cnt >= 1  && isaltsplice(motif)) 
         # psi = 0.99 && mle_ci( psi, inc_cnt, z=1.64 )
      else
         # NA is ignored.
      end
   else # do EM
      rec_spliced_em!(  )
   end
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


