
typealias NodeInt UInt32

immutable SGNodeTup
   gene::NodeInt
   node::NodeInt
end

typealias SGNodeSet Vector{SGNodeTup}

immutable Edges{K}
   left::Vector{SGNodeSet}
   right::Vector{SGNodeSet}
end

Base.convert{K}( ::Type{Edges{K}}, graphs::Vector{SpliceGraph} ) = build_edges( graphs, K )

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

function add_kmer_edge!{S <: NucleotideSequence}( kmers::Vector{SGNodeSet}, 
                                                  seq::S, l, r, entry::SGNodeTup )
   try
      println( "$(seq[(l-4):(l-1)]) + $(seq[l:r]) + $(seq[(r+1):(r+4)])" )
      curkmer = DNAKmer( seq[l:r] )
      ind = Int(UInt64(curkmer)) + 1
      if !isdefined(kmers, ind)
         kmers[ind] = SGNodeSet()
      end
      push!(kmers[ind], entry)
      kmers[ind] = sort(kmers[ind], by=x->x.gene) |> unique
   end
end

function build_edges( graphs::Vector{SpliceGraph}, k::Integer )
   @assert(1 <= k <= 32, "$k is out of bounds for Edges{K} set!")
   left  = Vector{SGNodeSet}(4^k)
   right = Vector{SGNodeSet}(4^k)

   for (i,g) in enumerate(graphs)
      for (j,n) in enumerate(g.nodeoffset)

         if is_edge( g.edgetype[j], true ) # left edge
            add_kmer_edge!( left, g.seq, n-k-2, n-3, SGNodeTup(i,j) )
         end
         if is_edge( g.edgetype[j], false ) # right edge
            add_kmer_edge!( right, g.seq, n,  n+k-1, SGNodeTup(i,j) )
         end

      end
   end
   Edges{k}(left, right)
end

