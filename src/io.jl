
#=
#
#                     Basic IO Functions...
#
=#

plug_write( io::BufOut, str::S; plug::Char='\t' ) where {S <: AbstractString} = (write( io, str ); write( io, plug  ))
plug_write( io::BufOut, str::Char; plug::Char='\t' ) = (write( io, str ); write( io, plug  ))

tab_write( io::BufOut, str::S ) where {S <: AbstractString} = plug_write( io, str, plug='\t' )
tab_write( io::BufOut, str::Char ) = plug_write( io, str, plug='\t' )

end_write( io::BufOut, str::S ) where {S <: AbstractString} = plug_write( io, str, plug='\n' )
end_write( io::BufOut, str::Char ) = plug_write( io, str, plug='\n' )

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
   write( io, COMPLEX_CHAR )
   write( io, string(complex) )
   tab && write( io, '\t' )
end

function conf_int_write( io::BufOut, conf_int::Tuple; tab=false, width=false, sig=4 )
   lo,hi = conf_int #unpack
   if width
      write( io, string( round(hi-lo, digits=sig) ) )
      write( io, '\t' )
   end
   write( io, string(lo) )
   write( io, ',' )
   write( io, string(hi) )
   tab && write( io, '\t' )
end

function edges_write( io::BufOut, edges::IntervalMap{ExonInt,C}, range::UnitRange ) where C <: ReadCounter
   first = true
   for i in edges
      if i.first in range && i.last in range
         !first && write( io, ',' )
         write( io, string(i.first) )
         write( io, '-' )
         write( io, string(i.last) )
         write( io, ':' )
         write( io, string(round(get(i.value), digits=1)) )
         first = false
      end
   end
end

#=
#
#                     SAM Format IO...
#
=#

function seq_write( io::BufOut, read::FASTQRecord; tab=false )
   for c in read.sequence
      write( io, convert(Char, c) )
   end
   tab && write( io, '\t' )
end

function qual_write( io::BufOut, read::FASTQRecord; qualoffset=33, tab=false )
   for i in read.quality
      write( io, convert(Char, i+qualoffset) )
   end
   tab && write( io, '\t' )
end

function sam_flag( align::SGAlignment, lib::GraphLib, ind, paired, first, is_pair_rc, is_supplemental )
   flag = UInt16(0)
   if paired
      flag |= 0x01
      flag |= 0x02
      if first
         flag |= 0x40
         lib.info[ ind ].strand || (flag |= 0x10)
      else
         flag |= 0x80
         lib.info[ ind ].strand || (flag |= 0x10)
         #is_pair_rc   && (flag ⊻= 0x20)
      end
   else
      lib.info[ ind ].strand || (flag ⊻= 0x10)
   end
   is_supplemental && (flag |= 0x100)
   flag
end

@inline function soft_pad( align::SGAlignment )
   pad = align.offsetread - 1
   pad > 0 ? string(pad) * "S" : ""
end

function cigar_string( align::SGAlignment, sg::SpliceGraph, strand::Bool, readlen=0 )
   curpos    = align.offset                                 # left most position
   running   = 0                                            # cur number of running matches
   total     = align.offsetread - 1                         # total positions accounted for
   cigar     = String[ soft_pad(align) ]                    # build cigar string here

   for idx in 1:length(align.path)

      matches = align.path[idx].score.matches + align.path[idx].score.mismatches
      total  += matches
      curpos += matches

      if idx < length( align.path )
         curi = align.path[idx].node
         nexti = align.path[idx+1].node
         intron = strand ? Int(sg.nodecoord[nexti]) - Int(sg.nodecoord[curi] +  sg.nodelen[curi] - 1) - 1 :
                                 Int(sg.nodecoord[curi])  - Int(sg.nodecoord[nexti] + sg.nodelen[nexti] - 1) - 1
         if intron >= 1
            push!( cigar, string(running + matches) * "M" )
            push!( cigar, string( intron ) * "N" )
            running = 0
            curpos  = sg.nodeoffset[nexti]
         else
            running += matches
         end
      else
         push!( cigar, string(running + matches) * "M" )
      end
   end
   if total < readlen
      soft = readlen - total
      push!( cigar, string( soft ) * "S" )
   end

   !strand && reverse!( cigar )
   join( cigar, "" ), curpos
end


# Write single SAM entry for one SGAlignment
function write_sam( io::BufOut, read::FASTQRecord, align::SGAlignment, lib::GraphLib;
                    mapq::Int=0, paired::Bool=false, fwd_mate::Bool=true, is_pair_rc::Bool=true, qualoffset::Int=33,
                    supplemental::Bool=false, tagstr::String="" )
   (align.isvalid && length(align.path) >= 1) || return
   geneind = align.path[1].gene
   strand  = lib.info[geneind].strand
   nodeind = strand ? align.path[1].node : align.path[end].node
   (align.path[end].node < nodeind || align.path[end].node < align.path[1].node) && return # TODO: allow circular SAM output
   sg = lib.graphs[geneind]
   cigar,endpos = cigar_string( align, sg, strand, length(read.sequence) )
   offset = strand ? (align.offset - sg.nodeoffset[nodeind]) : (sg.nodelen[nodeind] - (endpos - sg.nodeoffset[nodeind]))
   if !strand && !supplemental
      reverse_complement!(read)
   end
   tab_write( io, string(read) )
   tab_write( io, string( sam_flag(align, lib, geneind, paired, fwd_mate, is_pair_rc, supplemental) ) )
   tab_write( io, lib.info[geneind].name )
   tab_write( io, string( sg.nodecoord[nodeind] + offset ) )
   tab_write( io, string(mapq) )
   tab_write( io, cigar )
   tab_write( io, '*' )
   tab_write( io, '0' )
   tab_write( io, '0' )
   if !supplemental
      seq_write( io, read, tab=true )
      qual_write( io, read, qualoffset=qualoffset )
   else
      tab_write( io, '*' )
      write( io, '*' )
   end
   if tagstr != ""
      write( io, '\t' )
      write( io, tagstr )
   end
   write( io, '\n' )
end

# Efficient version for write_sam with multiple supplementary alignments
function indmax_score( vec::Vector{SGAlignment} )
   ind = 1
   max = score(vec[1])
   for i in 2:length(vec)
      cur = score(vec[i])
      if cur > max
         max = cur
         ind = i
      end
   end
   ind
end

function indmax_score( fwd::Vector{SGAlignment}, rev::Vector{SGAlignment} )
   ind = 1
   max = score(fwd[1], rev[1])
   for i in 2:length(fwd)
      cur = score(fwd[i], rev[i])
      if cur > max
         max = cur
         ind = i
      end
   end
   ind
end

# Write multiple SAM entries for a Vector of SGAlignment, highest score gets regular entry
# others get supplementary entries.
function write_sam( io::BufOut, read::FASTQRecord, alignvec::Vector{SGAlignment}, lib::GraphLib;
                    paired::Bool=false, fwd_mate::Bool=true, is_pair_rc::Bool=true, qualoffset::Int=33 )

   best = indmax_score( alignvec )
   mapq_val = 10 * length(alignvec)
   write_sam( io, read, alignvec[best], lib, mapq=mapq_val, paired=paired, fwd_mate=fwd_mate,
                        is_pair_rc=is_pair_rc, qualoffset=qualoffset, supplemental=false )
   for i in 1:length(alignvec)
      if i != best
         write_sam( io, read, alignvec[i], lib, mapq=mapq_val, paired=paired, fwd_mate=fwd_mate,
                        is_pair_rc=is_pair_rc, qualoffset=qualoffset, supplemental=true )

      end
   end
end

# Write multiple SAM entries for two vectors of SGAlignments for paired-end reads
function write_sam( io::BufOut, fwd::FASTQRecord, rev::FASTQRecord,
                    fwd_vec::Vector{SGAlignment}, rev_vec::Vector{SGAlignment}, lib::GraphLib;
                    paired::Bool=true, fwd_mate::Bool=true, is_pair_rc::Bool=true, qualoffset::Int=33 )

   best = indmax_score( fwd_vec, rev_vec )
   mapq_val = 10 * length(fwd_vec)
   write_sam( io, fwd, fwd_vec[best], lib,  mapq=mapq_val, paired=paired, fwd_mate=fwd_mate, is_pair_rc=is_pair_rc,
                                      qualoffset=qualoffset, supplemental=false )
   write_sam( io, rev, rev_vec[best], lib, mapq=mapq_val, paired=paired, fwd_mate=!fwd_mate, is_pair_rc=is_pair_rc,
                                      qualoffset=qualoffset, supplemental=false )
   for i in 1:length(fwd_vec)
      if i != best
         write_sam( io, fwd, fwd_vec[i], lib, mapq=mapq_val, paired=paired, fwd_mate=fwd_mate, is_pair_rc=is_pair_rc,
                                         qualoffset=qualoffset, supplemental=true )
         write_sam( io, rev, rev_vec[i], lib, mapq=mapq_val, paired=paired, fwd_mate=fwd_mate, is_pair_rc=is_pair_rc,
                                         qualoffset=qualoffset, supplemental=true )
      end
   end
end


function write_sam_header( io::BufOut, lib::GraphLib, head::Bool=true )
   refs = Dict{String,Int}()
   for gind in 1:length(lib.graphs)
      name = lib.info[gind].name
      len  = max( lib.graphs[gind].nodecoord... )
      if !haskey( refs, name ) || refs[name] < len
         refs[name] = len + 10000
      end
   end
   head && write( io, "@HD\tVN:1.0\tSO:unsorted\n" )
   for k in keys(refs)
      write( io, "@SQ\tSN:" )
      tab_write( io, k )
      write( io, "LN:" * string( refs[k] ) )
      write( io, '\n' )
   end
end

#=
#
#                     Event IO...
#
=#

function output_utr( io::BufOut, psi::Vector{Float64}, pgraph::Nullable{PsiGraph},
                     total_reads::Float64, motif::EdgeMotif, sg::SpliceGraph, node::Int,
                     info::GeneMeta; empty=false )
   st = motif == TXST_MOTIF ? node : node - 1
   en = st + length(psi) - 1
   if en < st
      error("$psi, $total_reads, $motif, $node, $info, $empty")
   end
   i = 1
   for n in st:en
      tab_write( io, info[1] )
      tab_write( io, string(n) )
      coord_write( io, info[2], sg.nodecoord[n], sg.nodecoord[n]+sg.nodelen[n]-1, tab=true )
      tab_write( io, info[3] )
      tab_write( io, convert(String, motif) )
      if empty # had to add this flag since we iterate through the TS/TE nodes
         write( io, "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" )
         i += 1
         continue
      end
      tab_write( io, string(psi[i]) )
      if !isnull( pgraph )
         conf_int  = beta_posterior_ci( psi[i], total_reads, sig=3 )
         conf_int_write( io, conf_int, tab=true, width=true)
         tab_write( io, string( round(total_reads, sigdigits=3) ) )
         complex_write( io, complexity( pgraph.value ), tab=true )
         tab_write( io, string( round( shannon_index( pgraph.value ), digits=4 ) ) )
         count_write( io, get(pgraph).nodes[i], get(pgraph).psi[i], tab=true )
         write( io, "NA\t" )
      else
         tab_write( io, "NA\tNA\tNA\tNA\tNA\tNA\tNA\t" )
      end
      write( io, "NA" )
      write( io, '\n' )
      i += 1
   end
   en
end

function output_psi( io::BufOut, psi::Float64, inc::Nullable{PsiGraph}, exc::Nullable{PsiGraph},
                     total_reads::Float64, conf_int::Tuple, motif::EdgeMotif,
                     sg::SpliceGraph, edges::IntervalMap{ExonInt,C},
                     node::Int, info::GeneMeta, bias ) where C <: ReadCounter

   sg.nodelen[node] == 0 && return
   tab_write( io, info[1] ) # gene
   tab_write( io, string(node) )
   coord_write( io, info[2], sg.nodecoord[node], sg.nodecoord[node]+sg.nodelen[node]-1, tab=true ) #coord
   tab_write( io, info[3] )
   tab_write( io, convert(String, motif) )
   tab_write( io, string(psi) )

   conf_int_write( io, conf_int, tab=true, width=true )
   tab_write( io, string( round(total_reads, sigdigits=3) ) )

   if !isnull( inc ) && !isnull( exc )
      complex_write( io, complexity( inc.value, exc.value ), tab=true )
      tab_write( io, string( round( shannon_index( inc.value, exc.value ), digits=4 ) ) )
   else
      write( io, "K0\t0.0\t" )
   end

   min,max = Inf,0
   if !isnull( inc )
      count_write( io, get(inc), psi, tab=true )
      min,max = get(inc).min,get(inc).max
   else
      write( io, "NA\t" )
   end

   if !isnull( exc )
      count_write( io, get(exc), 1.0-psi, tab=true )
      min = get(exc).min < min ? get(exc).min : min
      max = get(exc).max > max ? get(exc).max : max
   else
      write( io, "NA\t" )
   end

   edges_write( io, edges, min:max )
   write( io, '\n' )
end

function output_circular( io::BufOut, sg::SpliceGraph, sgquant::SpliceGraphQuant, info::GeneMeta )
   for (st,en) in keys(sgquant.circ)
      back_len = 0.0
      back_cnt = 0.0
      fore_cnt = get(sgquant.circ[(st,en)])
      for edg in intersect( sgquant.edge, (st, st) )
         if edg.first == st
            back_len += 1
            back_cnt += get(edg.value)
         end
      end
      total_reads = (fore_cnt + back_cnt)
      psi = fore_cnt / total_reads
      if !isnan(psi) && total_reads > 0
         tab_write( io, info[1] )
         tab_write( io, string(st) * "-" * string(en) )
         coord_write( io, info[2], sg.nodecoord[st]+sg.nodelen[st]-1, sg.nodecoord[en], tab=true )
         tab_write( io, info[3] )
         tab_write( io, "BS" )
         tab_write( io, string(psi) )
         conf_int  = beta_posterior_ci( psi, total_reads, sig=3 )
         conf_int_write( io, conf_int, tab=true, width=true)
         tab_write( io, string( round(total_reads, sigdigits=3) ) )
         write( io, "NA\tNA\tNA\tNA\tNA\n" )
      end
   end
end

function output_empty( io::BufOut, motif::EdgeMotif, sg::SpliceGraph, node::Int, info::GeneMeta )
   sg.nodelen[node] == 0 && return
   tab_write( io, info[1] )
   tab_write( io, string(node) )
   coord_write( io, info[2], sg.nodecoord[node], sg.nodecoord[node]+sg.nodelen[node]-1, tab=true )
   tab_write( io, info[3] )
   tab_write( io, convert(String, motif) )
   write( io, "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" )
end

function Base.string( is::BitSet )
   str = ""
   for n in is
      str *= string(n) * "-"
   end
   str[1:end-1]
end

function count_write( io::BufOut, nodestr, psi::Float64; tab=false )
   write( io, string(nodestr) )
   write( io, ':' )
   write( io, string(round(psi, digits=6)) )
   tab && write( io, '\t' )
end

function count_write( io::BufOut, nodestr, psi::String; tab=false )
   write( io, string(nodestr) )
   write( io, ':' )
   write( io, psi )
   tab && write( io, '\t' )
end

function count_write( io::BufOut, pgraph::PsiGraph, psi::Float64; tab=false )
   start = 1
   if psi == 1.0 && length(pgraph.nodes) > 1
      for i in 1:length(pgraph.nodes)-1
         write( io, string(pgraph.nodes[i]) )
         write( io, ',' )
      end
      start = length(pgraph.nodes)
   end
   for i in start:length(pgraph.nodes)
      curpsi = i > length(pgraph.psi) ? psi : pgraph.psi[i]
      count_write( io, pgraph.nodes[i], curpsi )
      (i < length(pgraph.nodes)) && write( io, "," )
   end
   tab && write( io, '\t' )
end

function output_psi_header( io::BufOut )
   tab_write( io, "Gene\tNode\tCoord\tStrand" )
   tab_write( io, "Type\tPsi\tCI_Width\tCI_Lo,Hi" )
   write( io, "Total_Reads\tComplexity\tEntropy\tInc_Paths\tExc_Paths\tEdges\n" )
end

function output_tpm_header( io::BufOut, idstr::String="Gene" )
   write( io, idstr )
   write( io, "\tTpM\tRead_Counts\n" )
end

function output_diff_header( io::BufOut )
   tab_write( io, "Gene\tNode\tCoord\tStrand" )
   write( io, "Type\tPsi_A\tPsi_B\tDeltaPsi\tProbability\tComplexity\tEntropy\n" )
end

function output_diff( io::BufOut, event, complex::Int, entropy::Float64,
                      psi_a::Float64, psi_b::Float64, deltapsi::Float64, prob::Float64, sig=5 )
   for i in 1:length(event)
      tab_write( io, event[i] )
   end
   tab_write( io, string(round(psi_a, sigdigits=sig)) )
   tab_write( io, string(round(psi_b, sigdigits=sig)) )
   tab_write( io, string(round(deltapsi, sigdigits=sig)) )
   tab_write( io, string(round(prob, sigdigits=sig)) )
   write( io, COMPLEX_CHAR )
   tab_write( io, string(complex) )
   tab_write( io, string(entropy) )
   write( io, '\n' )
end

#=
#
#                     Mapping Stats IO...
#
=#

function count_annotated_edges( sg::SpliceGraph, sgquant::SpliceGraphQuant )
   annoreads  = 0.0
   annocnt    = 0
   totalreads = 0.0
   totalcnt   = 0
   for edg in sgquant.edge
      totalreads += get(edg.value)
      totalcnt   += 1
      if edg in sg.annopath
         annoreads += get(edg.value)
         annocnt   += 1
      end
   end
   annoreads,annocnt,totalreads,totalcnt
end

function count_read_types( lib::GraphLib, graphq::GraphLibQuant )
   bodyreads = 0.0
   juncreads = 0.0
   annoreads = 0.0
   junccnt   = 0
   annocnt   = 0
   for i in 1:length(lib.graphs)
      bodyreads += sum(map(get, graphq.quant[i].node))
      areads,acnt,treads,tcnt = count_annotated_edges( lib.graphs[i], graphq.quant[i] )
      annoreads += areads
      juncreads += treads
      annocnt   += acnt
      junccnt   += tcnt
   end
   Int(round(bodyreads)),Int(round(juncreads)),Int(round(annoreads)),junccnt,annocnt
end

function output_stats( filename::String, lib::GraphLib, graphq::GraphLibQuant, param::AlignParam,
                       index::String, total::Int, mapped::Int, multi::Int, readlen::Int, ver::String )
   io = open( filename, "w" )
   stream = ZlibDeflateOutputStream( io )

   output_stats( stream, lib, graphq, param, index, total, mapped, multi, readlen, ver )

   close(stream)
   close(io)
end

function output_stats( io::BufOut, lib::GraphLib, graphq::GraphLibQuant, param::AlignParam,
                       index::String, total::Int, mapped::Int, multi::Int, readlen::Int, ver::String )
   write( io, "# -------- Whippet $ver ---------\n" )
   write( io, "Mapped_Index\t$index\n" )
   write( io, "Read_Length\t$readlen\n" )
   write( io, "Total_Reads\t$total\n" )
   write( io, "# -------- Mapping Stats ---------\n" )
   write( io, "Mapped_Reads\t$mapped\n" )
   write( io, "Mapped_Percent\t$(round((mapped/total)*100,sigdigits=2))%\n" )
   write( io, "# -------- Multi-loci Mapping ---------\n" )
   write( io, "Multimap_Reads\t$multi\n" )
   write( io, "Multimap_Percent\t$(round((multi/mapped)*100,sigdigits=2))%\n" )
   write( io, "# -------- Unspliced Reads ---------\n" )
   body,junc,anno,jcnt,acnt = count_read_types( lib, graphq )
   write( io, "Exon_Body_Reads\t$body\n" )
   write( io, "Exon_Body_Percent\t$(round((body/(body+junc))*100,sigdigits=2))%\n" )
   write( io, "# -------- Spliced Reads ---------\n" )
   write( io, "Junction_Reads\t$junc\n" )
   write( io, "Junction_Percent\t$(round((junc/(body+junc))*100,sigdigits=2))%\n" )
   write( io, "Junction_Number\t$jcnt\n" )
   write( io, "Annotated_Junc_Reads\t$anno\n" )
   write( io, "Annotated_Junc_Percent\t$(round((anno/junc)*100,sigdigits=2))%\n" )
   write( io, "Annotated_Junc_Number\t$acnt\n" )
   write( io, "Novel_Junc_Number\t$(jcnt-acnt)\n" )
   write( io, "Novel_Junc_Percent\t$(round(((jcnt-acnt)/jcnt)*100,sigdigits=2))%\n" )
end


function output_junctions( filename::String, lib::GraphLib, graphq::GraphLibQuant )
   io = open( filename, "w" )
   stream = ZlibDeflateOutputStream( io )

   output_junctions( stream, lib, graphq )

   close(stream)
   close(io)
end

function output_junctions( io::BufOut, lib::GraphLib, graphq::GraphLibQuant )
   function write_junctions( sg::SpliceGraph, sgquant::SpliceGraphQuant, i::Int )
      for edg in sgquant.edge
         if sg.edgetype[edg.first+1] in (EDGETYPE_LS, EDGETYPE_LL, EDGETYPE_LR) &&
            sg.edgetype[edg.last]    in (EDGETYPE_SR, EDGETYPE_RR, EDGETYPE_LR)

            if inall( edg, sg.annopath )
               write_junction( sg, edg, i, "CON_ANNO" )
            elseif edg in sg.annopath
               write_junction( sg, edg, i, "ALT_ANNO" )
            else # unique edge
               write_junction( sg, edg, i, "ALT_UNIQ" )
            end
         end
      end
   end

   function write_junction( sg::SpliceGraph, edg::IntervalValue, i::Int, str::String )
      tab_write( io, lib.info[i].name )
      if lib.info[i].strand
         donor    = sg.nodecoord[edg.first]+sg.nodelen[edg.first]-1
         acceptor = sg.nodecoord[edg.last]
      else
         donor    = sg.nodecoord[edg.last]+sg.nodelen[edg.last]-1
         acceptor = sg.nodecoord[edg.first]
      end
      tab_write( io, string(donor) )
      tab_write( io, string(acceptor) )
      plug_write( io, lib.info[i].gene, plug=':' )
      plug_write( io, string(edg.first), plug='-' )
      plug_write( io, string(edg.last), plug=':' )
      tab_write( io, str )
      tab_write( io, string(get(edg.value)) )
      end_write( io, lib.info[i].strand ? '+' : '-' )
   end

   write( io, "CHROM\tDONOR_COORD\tACCEPTOR_COORD\tGENEID:NODES:TYPE\tREADCOUNT\tSTRAND\n")
   for i in 1:length(lib.graphs)
      write_junctions( lib.graphs[i], graphq.quant[i], i )
   end
end


function output_exons( nodecoord::Vector{CoordInt},
                       nodelen::Vector{CoordInt},
                       tx::RefTx, strand::Bool )
   path = BitSet()
   # this may be `poor form`, but 256 is too big for default!
   # resize!(path.bits, 64) # Deprecated:  = zeros(UInt32,64>>>5)
   for i in 1:length(tx.acc)
      ind = collect(searchsorted( nodecoord, tx.acc[i], rev=!strand ))[end]
      push!( path, ind )
      cur = ind + (strand ? 1 : -1)
      while 1 <= cur <= length(nodecoord) &&
            nodecoord[cur] <= tx.don[i]
         push!( path, cur )
         cur = cur + (strand ? 1 : -1 )
      end
   end
   path
end

function Base.in( nodes::UnitRange, is::BitSet )
   for i in nodes
      if !(i in is)
         return false
      end
   end
   return true
end

function Base.in( nodes::UnitRange, v::Vector{BitSet} )
   for is in v
      if nodes in is
         return true
      end
   end
   return false
end

function output_exons( io::BufOut, sg::SpliceGraph, info::GeneInfo )
   for i in 1:length(sg.nodecoord)
      if isexonstart(sg.edgetype[i])
         for j in i:length(sg.nodecoord)
            if isexonstop(sg.edgetype[j+1])
               # print exon
               anno = i:j in sg.annopath ? "Y" : "N"
               tab_write( io, info.gene )
               if info.strand
                  coord_write( io, info.name, sg.nodecoord[i], sg.nodecoord[j] + sg.nodelen[j] - 1, info.strand ? "+" : "-", tab=true )
               else
                  coord_write( io, info.name, sg.nodecoord[j], sg.nodecoord[i] + sg.nodelen[i] - 1, info.strand ? "+" : "-", tab=true )
               end
               tab_write( io, join(i:j, ",") )
               end_write( io, anno )
            end
            if ishardedge(sg.edgetype[j+1])
               break
            end
         end
      end
   end
end

function output_exons( io::BufOut, lib::GraphLib )
   for i in 1:length(lib.graphs)
      output_exons( io, lib.graphs[i], lib.info[i] )
   end
end

function output_exons( filename::String, lib::GraphLib )
   io = open( filename, "w" )
   stream = ZlibDeflateOutputStream( io )

   write( stream, "Gene\tPotential_Exon\tWhippet_Nodes\tIs_Annotated\n" )
   output_exons( stream, lib )

   close(stream)
   close(io)
end
