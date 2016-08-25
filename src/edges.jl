
typealias NodeInt UInt32

immutable SGNode
   gene::NodeInt
   node::NodeInt
end

typealias SGNodeSet Vector{SGNode}

immutable Edges{K}
   left::Vector{SGNodeSet}
   right::Vector{SGNodeSet}
end

Base.(:<)( a::SGNode, b::SGNode ) = <( a.gene, b.gene )
Base.(:>)( a::SGNode, b::SGNode ) = >( a.gene, b.gene )
Base.(:(<=))( a::SGNode, b::SGNode ) = <=( a.gene, b.gene )
Base.(:(>=))( a::SGNode, b::SGNode ) = >=( a.gene, b.gene )

sortlt( a::SGNode, b::SGNode ) = a.gene == b.gene ? <( a.node, b.node ) : <( a.gene, b.gene )

Base.convert{K <: Integer}( ::Type{Edges{K}}, graphs::Vector{SpliceGraph} ) = build_edges( graphs, K )

# Is the EdgeType a connecting edge between two nodes?
# Check a donor splice site by:   is_edge( edge, true )
# and an acceptor splice site by: is_edge( edge, false )
# Return: Bool
function is_edge( edge::EdgeType, left::Bool )
   eint = convert(UInt8, edge)
   eint == 0x04 && return true
   left || (eint $= 0b011)
   (eint == 0x02 || eint == 0x05) && return true
   false
end

Base.intersect( arrA::Vector{SGNode}, arrB::Vector{SGNode} ) = intersect_sorted( arrA, arrB )

# Intersect two sorted arrays of SGNode by gene entry
# Return the list of intersected arrB entries
function intersect_sorted{T}( arrA::Vector{T}, arrB::Vector{T}; right=true )
   res = Vector{T}()
   i,j = 1,1
   while i <= length(arrA) && j <= length(arrB)
      if arrA[i] > arrB[j]
         j += 1
      elseif arrA[i] < arrB[j]
         i += 1
      else
         push!( res, right ? arrB[j] : arrA[i] )
         i += 1
         j += 1
      end
   end
   res
end 

kmer_index{T,K}( kmer::Kmer{T,K} ) = Int(UInt64(kmer)) + 1
kmer_index{T,K}( kmer::Bio.Seq.Kmer{T,K} ) = Int(UInt64(kmer)) + 1

kmer_index( seq::BioSequence ) = _kmer_index( seq )
kmer_index( seq::SGSequence  ) = _kmer_index( seq )

# calculate kmer index directly
function _kmer_index( seq )
   x     = UInt64(0)
   for nt in seq
      ntint = convert(UInt8, nt)
      if ntint > 0x03
         return 0 
      else
         x = x << 2 | ntint
      end
   end
   Int(x) + 1
end

function add_kmer_edge!{S <: NucleotideSequence}( kmers::Vector{SGNodeSet}, 
                                                  seq::S, l, r, left::Bool,
                                                  entry::SGNode )
   s = seq[l:r]
   ksize = r-l+1
   ind = 0 #default
   #println( "$(seq[(l-4):(l-1)]) + $(seq[l:r]) + $(seq[(r+1):(r+4)])" )
   try
      curkmer = SGKmer( s ) 
      ind = kmer_index(curkmer)
   catch
      abstr = AbstractString(s)
      ismatch( r"S|N", abstr ) && return(zero(UInt64))
      sub = replace( abstr, r"L|R", "" )
      #println("Caught $abstr replaced to $sub , $ksize")
      curl,curr = l-1,r+1
      while length(sub) < ksize && curl >= 1 && curr <= length(seq)
         sub = left ? AbstractString(seq[curl:curl]) * sub :
                      sub * AbstractString(seq[curr:curr])
         #println("Sub looks like $sub")
         sub = replace( sub, r"L|R", "" )
         curr += 1
         curl -= 1
      end
      try
        curkmer = SGKmer( s )
        ind = kmer_index(curkmer)
        #println( "Index size: $ind, kmer size: $(typeof(curkmer))" )
      end
   end
   if ind > 0
      if !isdefined(kmers, ind)
         kmers[ind] = SGNodeSet()
      end
      push!(kmers[ind], entry)
      kmers[ind] = sort( kmers[ind], lt=sortlt ) |> unique
   end
   UInt64(max(0, ind-1))
end

function build_edges( graphs::Vector{SpliceGraph}, k::Integer )
   @assert(1 <= k <= 32, "$k is out of bounds for Edges{K} set!")
   left  = Vector{SGNodeSet}(4^k)
   right = Vector{SGNodeSet}(4^k)

   for (i,g) in enumerate(graphs)
      for (j,n) in enumerate(g.nodeoffset)

         if is_edge( g.edgetype[j], true ) && isvalid( g.seq, (n-k-2):(n-3) ) # left edge
            lkmer = add_kmer_edge!( left, g.seq, n-k-2, n-3, true,  SGNode(i,j) )
            g.edgeleft[j] = SGKmer{k}(lkmer)
         end
         if is_edge( g.edgetype[j], false ) && isvalid( g.seq, n:(n+k-1) )# right edge
            rkmer = add_kmer_edge!( right, g.seq, n,  n+k-1, false, SGNode(i,j) )
            g.edgeright[j] = SGKmer{k}(rkmer)
         end

      end
   end
   Edges{k}(left, right)
end

