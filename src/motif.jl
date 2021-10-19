#using DataFrames
#using CSV
#using Glob
#using BioSequences
#using Statistics
#rbns = map(load_rbns, glob("*.tsv"))

encode_kmer( kmer::B ) where B <: BioSequence = Int(reinterpret(UInt64, DNAMer(kmer)))+1
encode_kmer( kmer::String ) = Int(reinterpret(UInt64, DNAMer(kmer)))+1
decode_kmer( index::Int, k::Int ) = String(reinterpret(DNAMer{k}, convert(UInt64, index - 1)))

encode_base( base::Char ) = Int(reinterpret(UInt64, DNAMer(base)))+1


mutable struct RNABindNSeq
   name::String
   conc::String
   kmer::Int
   gini::Float64
   data::Vector{Float64}
end

function gini_coefficient( arr::Vector{F} ) where F <: Number
   y    = !issorted(arr) ? sort(arr) : arr
   n    = length(y)
   wsum = 0.0
   cval = (n + 1) / n
   coef = (2 / n)
   for i in 1:length(y)
      wsum += i * y[i]
   end
   coef * wsum / sum(y) - cval
end

function load_rbns( tsv_file::String )
   rbp  = CSV.read(bufferedinput(tsv_file), DataFrame, delim='\t')
   gini = map(i->gini_coefficient(rbp[:,i]), 2:ncol(rbp))
   ind  = argmax(gini)
   name = replace(names(rbp)[1], r"\[|\]|_0nM" => "")
   leng = length(rbp[:,ind+1])
   klen = Int(floor(mean(map(x->length(x), rbp[:,1]))))
   data = zeros(4^klen)

   rbp[:,ind+1] = rbp[:,ind+1] .* (sortperm(rbp[:,ind+1]) / nrow(rbp))
   for i in 1:leng
   	kdex = encode_kmer(replace(rbp[i,1], "*" => ""))
   	data[kdex] = rbp[i,ind+1]
   end
   data = data / maximum(data)
   RNABindNSeq(name, 
   				names(rbp)[ind+1], 
   				klen, 
   				gini_coefficient(data), 
   				map(x->round(x, digits=3), data))
end

function score_cis( seq::B, rbp::RNABindNSeq ) where B <: BioSequence
	best = 0.0
   for i in 1:length(seq)-(rbp.kmer-1)
      ind  = encode_kmer( seq[i:i+rbp.kmer-1] )
      best = rbp.data[ind] > best ? rbp.data[ind] : best
   end
   round(best, digits=3)
end

function score_cis( seq::B, rbps::Vector{RNABindNSeq} ) where B <: BioSequence
   scores = zeros(length(rbps))
   for i in 1:length(scores)
      scores[i] = score_cis(seq, rbps[i])
   end
   scores
end

#=
	Rewritten in julia from MaxEntScan perl scripts --
	Yeo G and Burge C. Maximum entropy modeling of short sequence motifs with
	applications to RNA splicing signals. Journal of Computational Biology,
	2004, 11:377-94.
=#

const bgd_5   = [0.27, 0.23, 0.23, 0.27]
const cons1_5 = [0.004, 0.0032, 0.9896, 0.0032]
const cons2_5 = [0.0034, 0.0039, 0.0042, 0.9884]

const bgd_3   = [0.27, 0.23, 0.23, 0.27]
const cons1_3 = [0.9903, 0.0032, 0.0034, 0.0030]
const cons2_3 = [0.0027, 0.0037, 0.9905, 0.0030]

const hashd   = Dict('A' => '0', 'C' => '1', 'G' => '2', 'T' => '3')

#Deprecated, use encode_kmer instead, its faster
tr( seq ) = map(s -> get(hashd,s,s), seq)

function hashseq( seq )
	seq = tr(seq)
	hval = 1
	for (i,j) in enumerate(seq)
		hval += parse(Int, j) * 4 ^ (length(seq) - i)
	end
	return hval
end

load_matrix5( filename::String="../motif/maxent/score5_matrix.txt" ) = CSV.read(bufferedinput(filename), DataFrame, delim='\t', header=false)[:,2]
function load_matrix3( filename::String="../motif/maxent/score3_matrix.txt" )
	csv = CSV.read(bufferedinput(filename), DataFrame, delim='\t', header=false)
	mat = unstack(csv, :Column2, :Column3, allowduplicates=true)
   return mat[:,2:ncol(mat)]
end

#normalize_five( score::Float64, max_score=11.8 ) = normalize_score( score, max_score )
#normalize_three( score::Float64, max_score=15.52 ) = normalize_score( score, max_score )

function normalize_score( score, max_score=10.0, min_score=-3.0 )
   score = min( score, max_score )
   score = max( score, min_score )
   score += min_score * -1
   score /= max_score + min_score * -1
   round(score, digits=3)
end

score_five( seq::B, matrix=load_matrix5() ) where B <: BioSequence = score_five( string(seq), matrix )

function score_five( seq::String, matrix=load_matrix5() )
	#=
    Calculate 5' splice site strength
    (exon)XXX|XXXXXX(intron)
              **
    >>> round(score5('cagGTAAGT'), 2)
    10.86
    >>> round(score5('gagGTAAGT'), 2)
    11.08
    >>> round(score5('taaATAAGT'), 2)
    -0.12
    >>> matrix = load_matrix5()
    >>> round(score5('cagGTAAGT', matrix=matrix), 2)
    10.86
   =#
   if length(seq) != 9
      error("Wrong length of seq!")
   end

   # for key elements
	key = uppercase(seq[4:5])
	score = cons1_5[encode_base(key[1])] * cons2_5[encode_base(key[2])] / 
				(bgd_5[encode_base(key[1])] *   bgd_5[encode_base(key[2])])
	# for rest elements
	rest = uppercase(seq[1:3] * seq[6:9])
	rest_score = matrix[encode_kmer(rest)]
	# final score
	return log2(score * rest_score)
end

score_three( seq::B, matrix=load_matrix3() ) where B <: BioSequence = score_three( string(seq), matrix )

function score_three( seq, matrix=load_matrix3() )
	#=
	Calculate 3' splice site strength
	(intron)XXXXXXXXXXXXXXXXXXXX|XXX(exon)
	                          **
	>>> round(score3('ttccaaacgaacttttgtAGgga'), 2)
	2.89
	>>> round(score3('tgtctttttctgtgtggcAGtgg'), 2)
	8.19
	>>> round(score3('ttctctcttcagacttatAGcaa'), 2)
	-0.08
	>>> matrix = load_matrix3()
	>>> round(score3('ttccaaacgaacttttgtAGgga', matrix=matrix), 2)
	2.89
	=#
	# check length of fa
	if length(seq) != 23
	   error("Wrong length of seq!")
	end

	# for key elements
	key = uppercase(seq[19:20])
	score = cons1_3[encode_base(key[1])] * cons2_3[encode_base(key[2])] / 
				(bgd_3[encode_base(key[1])] *   bgd_3[encode_base(key[2])])

	# for rest elements
	rest = uppercase(seq[1:18] * seq[21:23])
	rest_score = 1.0
	rest_score *= matrix[1,hashseq(rest[1:7])]
	rest_score *= matrix[2,hashseq(rest[8:14])]
	rest_score *= matrix[3,hashseq(rest[15:length(rest)])]
	rest_score *= matrix[4,hashseq(rest[5:11])]
	rest_score *= matrix[5,hashseq(rest[12:18])]
	rest_score /= matrix[6,hashseq(rest[5:7])]
	rest_score /= matrix[7,hashseq(rest[8:11])]
	rest_score /= matrix[8,hashseq(rest[12:14])]
	rest_score /= matrix[9,hashseq(rest[15:18])]
	# final score
	return log2(score * rest_score)
end


