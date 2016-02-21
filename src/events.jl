
bitstype 8 EdgeMotif

Base.convert( ::Type{EdgeMotif}, motif::UInt8 ) = box(EdgeMotif, unbox(UInt8, motif ))
Base.convert( ::Type{UInt8}, motif::EdgeMotif ) = box(UInt8, unbox(EdgeMotif, motif ))

const NONE_MOTIF = convert(EdgeMotif, 0b1000)

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

const MOTIF_TABLE = fill(NONE_MOTIF, 128 )
      # Alt TxStart  SL SL
      MOTIF_TABLE[ 0b000000 ] = TXST_MOTIF
      # Alt PolyA    RS RS
      MOTIF_TABLE[ 0b011011 ] = TXEN_MOTIF

      # Alt FirstEx  SL LS
      MOTIF_TABLE[ 0b000010 ] = ALTF_MOTIF
      # Alt LastEx   SR RS
      MOTIF_TABLE[ 0b001011 ] = ALTL_MOTIF

      # RetainedInt  LL RR
      MOTIF_TABLE[ 0b101110 ] = RETI_MOTIF

      # Skipped      RR LL
      MOTIF_TABLE[ 0b110101 ] = SKIP_MOTIF
      #              RR LR
      MOTIF_TABLE[ 0b110100 ] = SKIP_MOTIF
      #              LR LL
      MOTIF_TABLE[ 0b100101 ] = SKIP_MOTIF
      #              LR LR
      MOTIF_TABLE[ 0b100100 ] = SKIP_MOTIF

      # Alt Donor    LL LL
      MOTIF_TABLE[ 0b101101 ] = ALTD_MOTIF
      #              LL LR
      MOTIF_TABLE[ 0b101100 ] = ALTD_MOTIF

      # Alt Acceptor RR RR
      MOTIF_TABLE[ 0b110110 ] = ALTA_MOTIF
      #              LR RR
      MOTIF_TABLE[ 0b100110 ] = ALTA_MOTIF

Base.convert(::Type{EdgeMotif}, tup::Tuple{EdgeType,EdgeType}) = MOTIF_TABLE[ (UInt8(tup[1]) << 3) | UInt8(tup[2]) ]

function process_sgquant( lib::GraphLib, graphq::GraphLibQuant )
   for g in 1:length(lib.graphs)
      bias_correct!( lib.graphs[g], graphq.quant[g] )
      calculate_psi( lib.graphs[g], graphq.quant[g] )
   end
end


function calculate_bias( sg::SpliceGraph, sgquant::SpliceGraphQuant )
   nodelen = relative_length(  )
end

function process_events( sg::SpliceGraph, sgquant::SpliceGraphQuant )
   # Heres the plan:
   # step through sets of edges, look for edge motifs, some are obligate calculations
   # others we only calculate psi if there is an alternative edge
   # Then calculate bias of the event 
   if length( sg.edges ) <= 2
      return 0
   end
   for i in 1:length(sg.edges)-1
      if obligate( sg.edges[i], sg.edges[i+1] )
            
   end
end
