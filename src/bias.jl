
abstract type ReadCounter end
abstract type BiasModel end

function pushto!(t::IntervalTrees.IntervalBTree{K, V, B}, key::I, value) where {K, V, B, I <: AbstractInterval{K}}
    return _push!(t.root, key, value)
end

function _push!(t::IntervalTrees.InternalNode{K, V, B}, key::AbstractInterval{K}, value) where {K, V, B}
    i = IntervalTrees.findidx(t, key)
    if 1 <= length(t) - 1 && key >= t.keys[i]
        return _push!(t.children[i+1], key, value)
    else
        return _push!(t.children[i], key, value)
    end
end

function _push!(t::IntervalTrees.LeafNode{K, V, B}, key::AbstractInterval{K}, value) where {K, V, B}
    i = IntervalTrees.findidx(t, key)
    if 1 <= i <= length(t) &&
        first(t.entries[i]) == first(key) && last(t.entries[i]) == last(key)
        push!( t.entries[i].value, value )
        return true
    end
    return false
end

function pushzero!( collection::Dict{K,V}, key, value ) where {K, V <: ReadCounter}
   if !haskey( collection, key )
      collection[key] = zero(V)
   end
   push!( collection[key], value )
end

function pushzero!( collection::IntervalMap{K,V}, key, value ) where {K <: Integer, V <: ReadCounter}
   if !haskey( collection, (key.first,key.last) )
      collection[key] = zero(V)
   end
   pushto!( collection, key, value )
end

function Base.get( cnt::R ) where R <: ReadCounter
   if !cnt.isadjusted
      error("Cannot get() read count from model thats not already adjusted!")
   else
      return cnt.count
   end
end

function Base.push!( cnt::R, value::F ) where {R <: ReadCounter, F <: AbstractFloat}
   if cnt.isadjusted
      cnt.count += value
   else
      cnt.count += value
      cnt.isadjusted = true
      #error("Cannot push! a float value to PrimerBiasCounter unless isadjusted is true!")
   end
   cnt
end

Base.identity( cnt::R ) where R <: ReadCounter = get( cnt )
Base.identity( cnt::IntervalValue{I,R} ) where R <: ReadCounter = get( cnt.value )

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

primer_adjust!( cnt::DefaultCounter, m::B ) where B <: BiasModel = (cnt.isadjusted = true)
gc_adjust!( cnt::DefaultCounter, m::B ) where B <: BiasModel = (cnt.isadjusted = true)

count!( mod::DefaultBiasMod, _ ) = 1.0
count!( mod::DefaultBiasMod, x, y ) = 1.0

# 5' Sequence bias model
# Implementation of Hansen et al, NAR 2010 

# hash by identity for speed
Base.hash( v::UInt16 ) = convert(UInt64, v)

function add_to_back!( back::Vector{Float64}, exp::ExpectedGC, weight::Float64 )
   for i in 1:length(exp)
      back[i] += exp[i] * weight
   end
   mod
end

primer_normalize!( mod::DefaultBiasMod )   = 1.0
gc_normalize!( mod::DefaultBiasMod, x, y ) = 1.0

struct JointBiasTemp
   primer::Vector{UInt16}
   gc::Vector{Int32}
end

struct JointBiasMod <: BiasModel
   fore::Vector{Float64}
   back::Vector{Float64}
   size::Int64
   backoffset::Int64
   backnpos::Int64
   foreoffset::Int64
   forenpos::Int64

   gcfore::Vector{Float64}
   gcback::Vector{Float64}
   gcsize::Int64
   gcwidth::Float64
   gcoffset::Int64
   gcincrement::Int64

   temp::JointBiasTemp

   function JointBiasMod( ksize::Int=6, backoff::Int=23, backn::Int=10, foreoff::Int=1, foren::Int=2,
                          gcoffset::Int=5, gcincrement::Int=10; gc_param = GCBiasParams(50, 0.05) )
      fore = zeros( 4^ksize )
      back = zeros( 4^ksize )
      primer_temp = zeros(UInt16, foren)

      gcfore = zeros(Float64, Int(ceil(1.0 / gc_param.width)+1))
      gcback = zeros(Float64, Int(ceil(1.0 / gc_param.width)+1))
      gctemp = zeros(Int32, 100)

      temp = JointBiasTemp( primer_temp, gctemp )

      new( fore, back, ksize, backoff, backn, foreoff, foren,
           gcfore, gcback, gc_param.length, gc_param.width, gcoffset, gcincrement, 
           temp )
   end
end

function primer_normalize!( mod::JointBiasMod )
   foresum = sum(mod.fore)
   backsum = sum(mod.back)
   @inbounds for i in 1:length(mod.fore)
      @fastmath mod.fore[i] = mod.fore[i] / foresum
      @fastmath mod.back[i] = mod.back[i] / backsum
   end
end

function primer_count!( mod::JointBiasMod, seq::BioSequence{A} ) where A <: BioSequences.Alphabet
   # 5 prime bias
   idx = 1
   for i in mod.foreoffset:(mod.foreoffset+mod.forenpos-1)
      hept = kmer_index_trailing(UInt16, seq[i:(i+mod.size-1)])
      mod.fore[hept+1] += 1.0
      mod.temp.primer[idx] = hept
      idx += 1
   end
   for i in mod.backoffset:(mod.backoffset+mod.backnpos-1)
      mod.back[kmer_index_trailing(UInt16, seq[i:(i+mod.size-1)])+1] += 1.0
   end

   return mod.temp
end

function gc_count!( mod::JointBiasMod, seq::BioSequence{A} ) where A <: BioSequences.Alphabet
   # GC fragment bias
   idx = 1
   for i in mod.gcoffset:mod.gcincrement:(length(seq)-mod.gcsize)
      gc = gc_content(seq[i:(i+mod.gcsize-1)])
      bin = Int(div( gc, mod.gcwidth )+1)
      mod.gcfore[bin] += 1.0
      mod.temp.gc[idx] = bin
      idx += 1
   end

   return mod.temp
end

function count!( mod::JointBiasMod, seq::BioSequence{A} ) where A <: BioSequences.Alphabet
   # 5 prime bias
   primer_count!( mod, seq )
   # GC fragment bias
   fill!(mod.temp.gc, zero(Int32))
   gc_count!( mod, seq )

   return mod.temp
end

function count!( mod::JointBiasMod, fwd::BioSequence{A}, rev::BioSequence{A} ) where A <: BioSequences.Alphabet
   # 5 prime bias
   primer_count!( mod, fwd )
   # GC fragment bias
   fill!(mod.temp.gc, zero(Int32))
   gc_count!( mod, fwd )
   gc_count!( mod, rev )

   return mod.temp
end

@fastmath value!( mod::JointBiasMod, hept::UInt16 ) = mod.back[hept+1] / mod.fore[hept+1]
@fastmath value!( mod::JointBiasMod, bin::Int ) = mod.fore[bin] > 0.0 ? mod.back[bin] / mod.fore[bin] : 1.0


mutable struct JointBiasCounter <: ReadCounter
   count::Float64
   map::Dict{UInt16,Int32}
   gc::Vector{Int32}
   isadjusted::Bool

   function JointBiasCounter( gc_param = GCBiasParams(50, 0.05) )
      gc = zeros(Int32, Int(ceil(1.0 / gc_param.width)+1))
      new(0.0, Dict{UInt16,Int32}(), gc, false)
   end
end

Base.zero(::Type{JointBiasCounter}) = JointBiasCounter()

function Base.push!( cnt::JointBiasCounter, temp::JointBiasTemp )
   for i in 1:length(temp.primer)
      push!( cnt, temp.primer[i] )
   end
   for i in 1:length(temp.gc)
      if temp.gc[i] > 0
         push!( cnt, temp.gc[i] )
      else
         return cnt
      end
   end
   cnt
end

function Base.push!( cnt::JointBiasCounter, bin::Int32 )
   cnt.gc[bin] += one(Int32)
   cnt
end

function Base.push!( cnt::JointBiasCounter, hept::UInt16 )
   if haskey( cnt.map, hept )
      cnt.map[hept] += 1
   else
      cnt.map[hept] = 1
   end
   cnt
end

adjust!( cnt::JointBiasCounter, mod::JointBiasMod ) = !cnt.isadjusted ? primer_adjust!( cnt, mod ) : gc_adjust!( cnt, mod )

@inline function primer_adjust!( cnt::JointBiasCounter, mod::JointBiasMod )
   for (k,v) in cnt.map
      cnt.count += value!( mod, k ) * v
   end
   cnt.count /= mod.forenpos
   cnt.isadjusted = true
end

@inline function gc_adjust!( cnt::JointBiasCounter, mod::JointBiasMod )
   weight = 0.0
   for i in 1:length(cnt.gc)
      weight += value!( mod, i ) * cnt.gc[i]
   end
   if sum(cnt.gc) > 0
      weight /= sum(cnt.gc)
      cnt.count *= weight
   end
   cnt
end
