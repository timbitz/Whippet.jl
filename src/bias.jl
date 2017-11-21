
abstract type ReadCounter end
abstract type BiasModel end

#This may introduce type-instability
#Though it may be necessary...
#Base.get( num::Float64 ) = num 

function pushzero!( collection::Dict{K,V}, key::K, value ) where {K, V <: ReadCounter}
   if !haskey( collection, key )
      collection[key] = zero(ReadCount)
   end
   push!( collection[key], value )
end

function pushzero!( collection::IntervalMap{K,V}, key, value ) where {K <: Integer, V <: ReadCounter}
   if !haskey( collection, (key.first,key.last) )
      collection[key] = zero(ReadCount)
   end
   push!( get(collection, key, IntervalValue(0,0,DEFAULTCOUNTER_ZERO)).value, value )
end

function Base.get( cnt::R ) where R <: ReadCounter
   if !cnt.isadjusted
      error("Cannot get() read count from model thats not already adjusted!")
   else
      return cnt.count
   end
end

struct DefaultBiasMod <: BiasModel end

mutable struct DefaultCounter <: ReadCounter
   count::Float64
   isadjusted::Bool

   DefaultCounter() = new( 1.0, true )
   DefaultCounter( i::Float64 ) = new( i, true )
end

const DEFAULTCOUNTER_ZERO = DefaultCounter(0.0)
const DEFAULTCOUNTER_ONE  = DefaultCounter(1.0)

Base.zero( ::Type{DefaultCounter} ) = DefaultCounter(0.0) #DEFAULTCOUNTER_ZERO
Base.one( ::Type{DefaultCounter} )  = DefaultCounter(1.0) #DEFAULTCOUNTER_ONE
Base.push!( cnt::DefaultCounter, value::Float64 ) = (cnt.count += value)

adjust!( cnt::DefaultCounter, m::B ) where B <: BiasModel = (cnt.isadjusted = true)
count!( mod::DefaultBiasMod, anything ) = 1.0

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

function normalize!( mod::PrimerBiasMod )
   foresum = sum(mod.fore)
   backsum = sum(mod.back)
   @inbounds for i in 1:length(mod.fore)
      @fastmath mod.fore[i] = mod.fore[i] / foresum
      @fastmath mod.back[i] = mod.back[i] / backsum
   end
end

function count!( mod::PrimerBiasMod, seq::BioSequence{A} ) where A <: BioSequences.Alphabet
   hept = kmer_index_trailing(UInt16, seq[1:mod.size])
   mod.fore[hept+1] += 1.0
   for i in mod.offset:(mod.offset+mod.npos)
      mod.back[kmer_index_trailing(UInt16, seq[i:(i+mod.size-1)])+1] += 1.0
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


# GC-Content Bias Correction

struct GCBiasParams
   length::Int64           # length offset
   width::Float64          # gc bin size
end

const ExpectedGC = Vector{Float16}

function ExpectedGC( seq::BioSequence{A}; gc_param = GCBiasParams(50, 0.05) ) where A <: BioSequences.Alphabet
   const bins = zeros(Float16, Int(div(1.0, gc_param.width)+1))
   for i in 1:length(seq)-gc_param.length+1
      gc = Int(div(gc_content(seq[i:i+gc_param.length-1]), gc_param.width)+1)
      @inbounds bins[gc] += 1.0
   end
   # normalize columns
   binsum = sum(bins)
   binsum == 0.0 && (return bins)
   for i in 1:length(bins)
      @fastmath bins[i] /= binsum
   end
   return bins
end

struct GCBiasMod <: BiasModel
   fore::ExpectedGC
   back::ExpectedGC
end

mutable struct GCBiasCounter <: ReadCounter
   count::Float64
   gc::ExpectedGC
   isadjusted::Bool
end

value!( mod::GCBiasMod, bin::Int ) = mod.back[bin] / mod.fore[bin]

function adjust!( cnt::GCBiasCounter, mod::GCBiasMod )
   for i in 1:length(cnt.gc)
      cnt.count += value!( mod, i ) * cnt.gc[i]
   end
   cnt.isadjusted = true
end

function Base.push!( cnt::GCBiasCounter, gc::F ) where F <: AbstractFloat
   const bin = div( gc, 1.0 / length(cnt.gc) )
   gc[bin] += 1.0
end

#=
# Combined Bias Correction

mutable struct WhippetBiasCounter <: ReadCounter
   count::Float64
   map::Dict{UInt16,Int32}
   gc::ExpectedGC
   isadjusted::Bool   
end
=#

