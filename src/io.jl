
tab_write{S <: AbstractString}( io::BufOut, str::S ) = (write( io, str ); write( io, '\t'  ))
tab_write( io::BufOut, str::Char ) = (write( io, str ); write( io, '\t'  ))

function coord_write( io::BufOut, chr, first, last; tab=false )
   write( io, chr   )
   write( io, ':'   )
   write( io, string(first) )
   write( io, '-'   )
   write( io, string(last)  )
   tab && write( io, '\t' )
end

function coord_write( io::BufOut, chr, first, last, strand; tab=false )
   coord_write( io, chr, first, last )
   write( io, ':'    )
   write( io, strand )
   tab && write( io, '\t' )
end

function complex_write( io::BufOut, complex::Int; tab=false )
   write( io, 'C' )
   write( io, string(complex) )
   tab && write( io, '\t' )
end

function conf_int_write( io::BufOut, conf_int::Tuple; tab=false, width=false, sig=4 )
   lo,hi = conf_int #unpack
   if width
      write( io, string( signif(hi-lo, sig) ) )
      write( io, '\t' )
   end
   write( io, string(lo) )
   write( io, ',' )
   write( io, string(hi) )
   tab && write( io, '\t' )
end

function seq_write( io::BufOut, read::SeqRecord; tab=false )
   for c in read.seq
      write( io, convert(Char, c) )
   end
   tab && write( io, '\t' )
end

function qual_write( io::BufOut, read::SeqRecord, qualoffset=64; tab=false )
   for i in read.metadata.quality
      write( io, convert(Char, i+qualoffset) )
   end
   tab && write( io, '\t' )
end

function sam_flag( align::SGAlignment, lib::GraphLib, ind )
   flag = UInt16(0)
   lib.info[ ind ].strand || (flag |= 0x10)
   flag   
end


# Lets build a pseudo spliced CIGAR string from an SGAlignment
function cigar_string( align::SGAlignment, sg::SpliceGraph, readlen=align.matches )
   matchleft = align.matches
   cigar = ""
   curpos = align.offset
   leftover = 0
   total = 0
   # step through nodes in the path
   for idx in 1:length( align.path )
      const i = align.path[idx].node
      i <= length(sg.nodeoffset) || return cigar # this shouldn't happen
      # do the remaining alignment matches fit into the current node width
      if matchleft + curpos <= sg.nodeoffset[i] + sg.nodelen[i] 
         cigar *= string( min( readlen - total, matchleft + leftover ) ) * "M"
         total += matchleft + leftover
         matchleft = 0
         leftover = 0
      else # they don't fit -->
         curspace = (sg.nodeoffset[i] + sg.nodelen[i] - 1) - curpos
         matchleft -= curspace
         curpos += curspace
         # is the read spliced and is there another node left
         if idx < length( align.path )
            nexti = align.path[idx+1].node
            nexti <= length(sg.nodeoffset) || return cigar
            curpos = sg.nodeoffset[ nexti ]
            const intron = sg.nodecoord[nexti] - (sg.nodecoord[i] + sg.nodelen[i] - 1)
            if intron > 0
               cigar *= string( min( readlen-total, curspace ) ) * "M" * string( intron ) * "N"
               total += curspace
            else
               leftover += curspace
            end
         end
      end
   end
   
   if matchleft + leftover > 0
      cigar *= string( min( readlen - total, matchleft + leftover) ) * "M"
      total += matchleft + leftover
   end
   if total < readlen
      soft = readlen - total
      cigar *= string( soft ) * "S"
   end
   cigar
end

function write_sam( io::BufOut, read::SeqRecord, align::SGAlignment, lib::GraphLib; mapq=0 )
   const geneind = align.path[1].gene
   const nodeind = align.path[1].node
   const sg = lib.graphs[geneind] 
   tab_write( io, read.name )
   tab_write( io, string( sam_flag(align, lib, geneind) ) )
   tab_write( io, lib.info[geneind].name )
   tab_write( io, string( sg.nodecoord[nodeind] + (align.offset - sg.nodeoffset[nodeind]) ) )
   tab_write( io, string(mapq) )
   tab_write( io, cigar_string( align, sg, length(read.seq) ) )
   tab_write( io, '*' )
   tab_write( io, '0' )
   tab_write( io, '0' )
   seq_write( io, read, tab=true )
   qual_write( io, read )
   write( io, '\n' )
end

function write_sam_header( io::BufOut, lib::GraphLib )
   refs = Dict{ASCIIString,Int}()
   for gind in 1:length(lib.graphs)
      name = lib.info[gind].name
      len  = max( lib.graphs[gind].nodecoord... )
      if !haskey( refs, name ) || refs[name] < len
         refs[name] = len + 10000
      end
   end
   write( io, "@HD\tVN:1.0\tSO:unsorted\n" )
   for k in keys(refs)
      write( io, "@SQ\tSN:" )
      tab_write( io, k )
      write( io, "LN:" * string( refs[k] ) )
      write( io, '\n' )
   end
end

### EVENT PRINTING
function output_utr( io::BufOut, psi::Vector{Float64}, pgraph::Nullable{PsiGraph}, 
                     total_reads::Float64, motif::EdgeMotif, sg::SpliceGraph, node::Int,
                     info::GeneMeta; empty=false )
   st = motif == TXST_MOTIF ? node : node - 1
   en = st + length(psi) - 1
   i = 1
   for n in st:en
      tab_write( io, info[1] )
      coord_write( io, info[2], sg.nodecoord[n], sg.nodecoord[n]+sg.nodelen[n]-1, tab=true )
      tab_write( io, info[3] )
      tab_write( io, convert(ASCIIString, motif) )
      if empty # had to add this flag since we iterate through the TS/TE nodes
         write( io, "NA\tNA\tNA\tNA\tNA\tNA\tNA\n" )
         continue
      end
      if !isnull( pgraph )
         complex_write( io, complexity( pgraph.value ), tab=true )
      else
         tab_write( io, "NA" )
      end
      tab_write( io, string(psi[i]) )
      if !isnull( pgraph )
         conf_int  = beta_posterior_ci( psi[i], total_reads, sig=3 )
         conf_int_write( io, conf_int, tab=true, width=true)
         tab_write( io, string( signif(total_reads, 3) ) )
         count_write( io, get(pgraph).nodes[i], get(pgraph).count[i], get(pgraph).length[i], tab=true )
      else
         tab_write( io, "NA\tNA\tNA\tNA" )
      end 
      write( io, "NA" )
      write( io, '\n' )
      i += 1
   end
   en
end

function output_psi( io::BufOut, psi::Float64, inc::Nullable{PsiGraph}, exc::Nullable{PsiGraph},
                     total_reads::Float64, conf_int::Tuple, motif::EdgeMotif, sg::SpliceGraph, node::Int,
                     info::GeneMeta, bias )

   tab_write( io, info[1] ) # gene
   coord_write( io, info[2], sg.nodecoord[node], sg.nodecoord[node]+sg.nodelen[node]-1, tab=true ) #coord
   tab_write( io, info[3] )
   tab_write( io, convert(ASCIIString, motif) )
   if !isnull( inc ) && !isnull( exc )
      complex_write( io, complexity( inc.value, exc.value ), tab=true )
   else
          tab_write( io, "NA" )
   end
   tab_write( io, string(psi) )

   if !isnull( inc ) && !isnull( exc )
      conf_int_write( io, conf_int, tab=true, width=true )
      tab_write( io, string( signif(total_reads, 3) ) )
      count_write( io, get(inc), tab=true )
      count_write( io, get(exc) )
   else
      write( io, "NA\tNA\tNA\tNA\tNA" )
   end

   write( io, '\n' )
end

function output_circular( io::BufOut, sg::SpliceGraph, sgquant::SpliceGraphQuant, info::GeneMeta )
   for (st,en) in keys(sgquant.circ)
      back_len = 0.0
      back_cnt = 0.0
      fore_cnt = sgquant.circ[(st,en)]
      for edg in intersect( sgquant.edge, (st, st) )
         if edg.first == st
            back_len += 1
            back_cnt += edg.value
         end
      end
      total_reads = (fore_cnt + back_cnt)
      psi = fore_cnt / total_reads
      if !isnan(psi) && total_reads > 0
         tab_write( io, info[1] )
         coord_write( io, info[2], sg.nodecoord[st]+sg.nodelen[st]-1, sg.nodecoord[en], tab=true )
         tab_write( io, info[3] )
         tab_write( io, "BS" )
         tab_write( io, "C1" )
         tab_write( io, string(psi) )
         conf_int  = beta_posterior_ci( psi, total_reads, sig=3 )
         conf_int_write( io, conf_int, tab=true, width=true)
         tab_write( io, string( signif(total_reads, 3) ) )
         write( io, "NA\tNA\n" )
      end
   end
end

function output_empty( io::BufOut, motif::EdgeMotif, sg::SpliceGraph, node::Int, info::GeneMeta )
   tab_write( io, info[1] )
   coord_write( io, info[2], sg.nodecoord[node], sg.nodecoord[node]+sg.nodelen[node]-1, tab=true )
   tab_write( io, info[3] )
   tab_write( io, convert(ASCIIString, motif) )
   write( io, "NA\tNA\tNA\tNA\tNA\tNA\tNA\n" )
end

function count_write( io::BufOut, nodestr, countstr, lengstr; tab=false )
   write( io, string(nodestr) )
#=   write( io, "-" )
   write( io, string(countstr) )
   write( io, "(" )
   write( io, string(lengstr) )
   write( io, ")" ) =#
   tab && write( io, '\t' )
end

function count_write( io::BufOut, pgraph::PsiGraph; tab=false )
   for i in 1:length(pgraph.nodes)
      count_write( io, pgraph.nodes[i], pgraph.count[i], pgraph.length[i] )
      (i < length(pgraph.nodes)) && write( io, "," )
   end
   tab && write( io, '\t' )
end

function output_psi_header( io::BufOut )
   tab_write( io, "Gene\tNode\tStrand" )
   tab_write( io, "Type\tComplexity\tPsi\tCI_Width\tCI_Lo,Hi" )
   write( io, "Total_Reads\tInc_Paths\tExc_Paths\n" )
end

function output_tpm_header( io::BufOut )
   write( io, "Gene\tTpM\tRead_Counts\n" )
end

function output_diff_header( io::BufOut )
   tab_write( io, "Gene\tNode\tStrand" )
   write( io, "Type\tComplexity\tPsi_A\tPsi_B\tDeltaPsi\tProbability\n" )
end

function output_diff( io::BufOut, event, complex::Int, psi_a::Float64, psi_b::Float64, deltapsi::Float64, prob::Float64, sig=5 )
   for i in 1:length(event)
      tab_write( io, event[i] )
   end
   write( io, 'C' )
   tab_write( io, string(complex) )
   tab_write( io, string(signif(psi_a, sig)) )
   tab_write( io, string(signif(psi_b, sig)) )
   tab_write( io, string(signif(deltapsi, sig)) )
   write( io, string(signif(prob, sig)) )
   write( io, '\n' )
end
