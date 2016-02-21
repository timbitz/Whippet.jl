
# OBLIGATE MOTIFS (provided the cognate nodes both have reads):

# Alt TxStart                  SL .. SL .. +
#                 
# Alt Polyadenylation     + .. RS .. RS
#
# Alt First Exon         SL .. LS .. LR/LL 
#
# Alt Last Exon       RR/LR .. SR .. RS 


# REQUIRES SPANNING EDGE:

# Retained Intron              LL .. RR
#
# Skipped Exon              RR/LR .. LR/LL
#
# Alt 5' Splice Site           LL .. LL/LR
#
# Alt 3' Spilce Site        RR/LR .. RR
#

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
