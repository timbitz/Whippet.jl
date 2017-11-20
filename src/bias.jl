
abstract type ReadCounter end
abstract type BiasModel end

struct DefaultBiasMod <: BiasModel end

mutable struct DefaultCounter <: ReadCounter
   count::Float64
   isadjusted::Bool

   DefaultCounter() = new( 1.0, true )
end

adjust!( cnt::DefaultCounter, m::B ) where B <: BiasModel = (cnt.isadjusted = true)
Base.push!( cnt::DefaultCounter, value::Float64 ) = (cnt.count += value)


# 5' Sequence bias model
# Implementation of Hansen et al, NAR 2010 

struct PrimerBiasMod <: BiasModel
   fore::Vector{Float64}
   back::Vector{Float64}
   size::Int64
   offset::Int64
   npos::Int64

   function PrimerBiasMod( ksize::Int=7, npos::Int=6, offset::Int=24 )
      fore = zeros( 4^ksize )
      back = zeros( 4^ksize )
      new( fore, back, ksize, offset, npos )
   end
end

function count!( mod::PrimerBiasMod, seq::BioSequence )
   hept = kmer_index_trailing(UInt16, seq[1:mod.size])
   mod.fore[hept+1] += 1.0
   for i in mod.offset:(mod.offset+mod.npos)
      mod.back[kmer_index_trailing(UInt16, seq[i:(i+mod.npos-1)])+1] += 1.0
   end
   return hept
end

@fastmath value!( mod::PrimerBiasMod, hept::UInt16 ) = mod.back[hept] / mod.fore[hept]

# hash by identity for speed
Base.hash( v::UInt16 ) = v

mutable struct PrimerBiasCounter <: ReadCounter
   count::Float64
   map::Dict{UInt16,Int32}
   isadjusted::Bool

   PrimerBiasCounter() = new( 0.0, Dict{UInt16,Int32}(), false )
end

@inline function adjust!( cnt::PrimerBiasCounter, mod::PrimerBiasMod )
   for (k,v) in cnt.map
      cnt.count += value!( mod, k ) * v
   end
   cnt.isadjusted = true
end

function Base.push!( cnt::PrimerBiasCounter, hept::UInt16 )
   if haskey( cnt.map, hept )
      cnt.map[hept] += 1
   else
      cnt.map[hept] = 1
   end
end


