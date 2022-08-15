#Run this script before everything else. Ignore the errors since its just the hashtags. Or just spam shift+enter
#this is just a few functions that slightly help with exploring the .fna file and getting kmer spectra (slow and useless)

#remember to have this vers of ngsu:  Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="Missing-LongCharSeq-fix", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))

#module SimpleExplore

#export seqLenV, avgSeqLen,recordCount,avgRecLen, seqMode, seqSd, flipDict, dictToVec,convertSeq, querySplit, kmerEucDist

using Plots, BioSequences,FASTX, NextGenSeqUtils, WebBlast, Distances #PyPlot messes with soem stuff sometimes
using CSV, DataFrames, Statistics, ProgressMeter, StringDistances, SeqUMAP, StatsBase, DelimitedFiles, StatsPlots

import Base.length

#This is also in ExactMatch.jl
function length(reader::FASTX.FASTA.Reader)
    len = 0
    for record in reader
        len += FASTX.FASTA.seqlen(record)
    end
    return len
end

#all nessecary genomic data:
VicPac = "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna"
fVicPac = FASTA.sequence((first(open(FASTA.Reader,VicPac))))
LA = "C:/Users/lu_41/Desktop/Sofo Prok/genbank_pull_lama_and_alpaca.fasta"
reader = open(FASTA.Reader,LA) #has 3 records
sg = dna""
for record in reader
    sg = sg*FASTA.sequence(record)
end
close(reader)
V3 = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/final/V3.fasta"
AlpacaV = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/AlpacaV.fasta"

#Functions to do with lengths and counts of sequences/assemblies: (work in progress)
function seqLenV(reader::FASTX.FASTA.Reader{})
    lengths = Int64[]
    for record in reader
        l = FASTX.FASTA.seqlen(record)
        push!(lengths,l)
    end
    return lengths
end

function seqLenV(path::String)
    seqLenV(open(FASTA.Reader, path))
end

function avgSeqLen(lengths::Vector{Int64})
    return sum(lengths)/length(lengths)
end

#counts the number of records in a reader.
function recordCount(reader::FASTX.FASTA.Reader)
    c = 0
    for record in reader
        c += 1
    end
    return c
end

function recordCount(path::String)
    reader = open(FASTA.Reader, path)
    c = 0
    for record in reader
        c += 1
    end
    return c
end


#BIG PROBLEM: I THINK DOING OPEN(FASTA.READER) DOESNT WORK IN FUNCTIONS.

function avgRecLen(reader::FASTX.FASTA.Reader)
    reader = open(FASTA.Reader, AlpacaV)
    Alen = 0
    for record in reader
        Alen += FASTX.FASTA.seqlen(record)
    end
    return Int64(round((Alen)/(recordCount(open(FASTA.Reader,AlpacaV)))))
end

function avgRecLen(path::String)
    reader = open(FASTA.Reader, path)
    Alen = 0
    for record in reader
        Alen += FASTX.FASTA.seqlen(record)
    end
    return (Alen)/recordCount(open(FASTA.Reader,AlpacaV))
end

function seqMode(lengths)
    return maximum(lengths) - minimum(lengths)
end

function seqSd(reader::FASTX.FASTA.Reader{})
    lengths = SeqLengths(reader)
    avg = sum(lengths)/length(lengths)
    all = []
    for length in lengths
        a = abs(length-avg)
        push!(all,a)
    end
    l = sum(all)
    return sum(all)/sqrt(length(lengths))
end


#Very important functions to do with Dictionary and vector conversions. ApproxMatch depends on these functions!

#function to flip all keys and values in a dictionary that will become useful later.
function flipDict(dict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    nd = Dict{Int64,LongSequence{DNAAlphabet{4}}}()
    for key in dict
        k = first(key)
        v = last(key)
        nd[v] = k
    end
    return nd
end

#kmer to interger and the reverse. Much more methods are possible (e.g. removing need to dictionary) but ill implement those later.
# need to do: dict vec, vec vec, vec dict methods.
function convertSeq(Counts::Dict{LongSequence{DNAAlphabet{4}},Int64}, KmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    ans = Dict{Int64,Int64}()
    for i in Counts
        i = first(i)
        ans[KmerDict[i]] = Counts[i]
    end
    return ans
end

function convertSeq(Counts::Dict{Int64,Int64}, KmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    dict = flipDict(KmerDict)
    ans = Dict{LongSequence{DNAAlphabet{4}},Int64}()
    for index in Counts #remember the key is the kmer and value is the count.
        kmerInt = first(index)
        ans[dict[kmerInt]] = last(index)
    end
    return ans
end

function convertSeq(IMGTref::Dict{LongSequence{DNAAlphabet{4}},Float64}, KmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    ans = Dict{Int64,Float64}()
    for i in IMGTref
        i = first(i)
        ans[KmerDict[i]] = IMGTref[i]
    end
    return ans
end


#function to convert kmerdict to vector while preserving order. More methods can be added later.
function dictToVec(KmerDict::Dict{Int64, Int64})
    ans = Int64[]
    k = length(first(first(KmerDict)))
    for i in 1:length(KmerDict)
        push!(ans, KmerDict[i])
    end
    return ans
end

function dictToVec(KmerDict::Dict{Int64, Float64})
    ans = Float64[]
    k = length(first(first(KmerDict)))
    for i in 1:length(KmerDict)
        push!(ans, KmerDict[i])
    end
    return ans
end

#In testing, ORder IS PRESERVED!

function dictToVec(KmerDict::Dict{LongSequence{DNAAlphabet{4}}, Float64})
    ans = Float64[]
    k = length(first(first(KmerDict)))
    for i in KmerDict
        push!(ans, KmerDict[first(i)])
    end
    return ans
end


#Splitting a sequence into constituent kmers in O(n) time:
#bens pcakage also had something like this idk if its faster
function querySplit(seq::LongSequence{DNAAlphabet{4}}, k::Int64, Dictionary::Bool = true)
    kmer_vector::Vector{LongSequence{DNAAlphabet{4}}} = [seq[i:i+k-1] for i in 1:length(seq)-k+1]
    if Dictionary == true
        kmer_dict = Dict{LongSequence{DNAAlphabet{4}}, Vector{Int64}}()
        for (i, kmer) in enumerate(kmer_vector)
            if haskey(kmer_dict, kmer)
                push!(kmer_dict[kmer], i)
            else
                kmer_dict[kmer] = [i]
            end
        end
        return kmer_dict
    elseif Dictionary == false
        return kmer_vector
    end
end


#euclidian distance. Will only be computed once in final function. ALthough now I uses the Distances.Ecuclidean function and its blazing fast. They probably use lots of bit and optimization tricks.
#I should test vector and dictionary versions for speed
#also i realized that a vector is completely doable but i wont implement it now.
function kmerEucDist(ref::Dict{Int64,Int64}, seq::Dict{Int64,Int64}) #assuming both are kmer freqs with same length.
    sqrdiffs = 0
        for key in ref
            diff = (last(key) - last(seq[first(key)]))^2
            sqrdiffs += diff
            dist = sqrt(sqrdiffs)
        end
        return dist
    end

#note that the distances package is faster

#old function from seconddraft might be useful just counts all kmers
#alternative that does count zeroes w kmerDict. note!!! its also possible to make a version that just generates a kmerdict full of zeroes!
function recordKCount(k::Int64, seq::LongSequence{DNAAlphabet{4}},kmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64}) #doesnt work for k=1 lol but that can be easily fixed later
    k = k-1
    rv = kmerDict
    for key in rv
        rv[first(key)] = 0
    end
    for i in 1:length(seq)-k
        subseq = seq[i:i+k]
        if get(rv,subseq,nothing) != nothing
            rv[subseq] += 1
        end
    end
    return rv
end

#fc = recordKCount(6,fVicPac,sixMerDict)

#version that does whole reader but not count zeroes. It takes way too long to run. already 5 minutes for k = 3
#there are actual good kmer counting algorithms but this is a naive implementation
function ReaderKCount(k::Int64, reader::FASTX.FASTA.Reader, kmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    rv = kmerDict
    for key in rv
        rv[first(key)] = 0
    end
    c = dna""
    if k > 1
        #k -= 1
        for record in reader
            seq = FASTA.sequence(record)
            Aseq = c*seq
            c = seq[end-k+2:end]
            for i in 1:length(Aseq)-k+1
                subseq = Aseq[i:i+k-1]
                if haskey(rv,subseq) #Fix by checking if its in the kmer dict.
                    rv[subseq] += 1
                end
            end
        end
    end
    return rv
end

#VP = open(FASTA.Reader,VicPac)
#@time Vp6k = ReaderKCount(6,VP,genKmers(6))

#maybe kmer spectra that fills a vector. KAT.jl has cooler ways to do this. Also this doesnt work lol
#problem: for something like VicPac theres too many different values. Perhaps making a histogram with each
function kmerSpectra(k::Int64, KCDict::Dict{LongSequence{DNAAlphabet{4}},Int64}, rm::Bool = false)
    ans = fill(0,(2^((4*k-1)))-(2^((2*k)-1)))
    for key in KCDict
        ans[last(key)+1] += 1
    end
    if rm
        Ransvec = reverse(ans)
        stop = false
        for i in Ransvec
            if stop == false
                if i == 0
                    pop!(ans)
                else
                    stop = true #there should be a better way
                end
            else
                break
            end
        end
    end
    Plots.scatter(ans,label=nothing)
    xlabel!(string(k)*"-mer frequency + 1")
    ylabel!("count")
end

#vkmerSpectra(6,Vp6k)
#kmerSpectra(6,fc)
#end

function oldKmerSpectra(k::Int64, KCDict::Dict{LongSequence{DNAAlphabet{4}},Int64},typ::String="bar")
    ans = Dict{Int64,Int64}()
    @showprogress "initializing frequency array..." for key in KCDict
        if haskey(ans,last(key))
            ans[last(key)] += 1
        else
            ans[last(key)] = 1
        end
    end
    ansvec = fill(0,(2^((4*k-1)))-(2^((2*k)-1))) #this is to minimize the amount of terms to put in the vector
    for key in ans
        ansvec[first(key)+1]=last(key)
    end
    Ransvec = reverse(ansvec)
    stop = false
    for i in Ransvec
        if stop == false
            if i == 0
                pop!(ansvec)
            else
                stop = true #there should be a better way
            end
        else
            break
        end
    end
    if typ == "bar"
        Plots.bar(ansvec,label=nothing)
    else
        Plots.scatter(ansvec,label=nothing)
    end
    xlabel!(string(k)*"-mer frequency + 1")
    ylabel!("count")
end

#oldKmerSpectra(6,fc)
#oldKmerSpectra(6,Vp6k)

##looking at FASTQ sequences with weird lengths(also can be done for fasta.):
function inspectSeqLen(reader::FASTX.FASTQ.Reader, lengths::UnitRange{Int64}, amount::Int64 = 10, A::Bool = false)
    c = 0
    if !A
        for record in reader
            if c <= amount
                len = FASTQ.seqlen(record)
                if len <= last(lengths) && len >= first(lengths)
                    println(record)
                    c += 1
                end
            end
        end
    else
        for record in reader
            if c <= amount
                len = FASTQ.seqlen(record)
                if len <= last(lengths) && len >= first(lengths)
                    println(">"*FASTQ.identifier(record)*" | length = "*string(len))
                    println(FASTQ.sequence(record))
                    c += 1
                end
            end
        end
    end
end

function howManySeq(reader::FASTX.FASTQ.Reader, n::Int64)
    lessthan = 0
    morethan = 0
    for record in reader
        len = FASTQ.seqlen(record)
        if len <= n
            lessthan += 1
        else
            morethan += 1
        end
    end
    return string(lessthan)*" below or equal to thr, "*string(morethan)*" above thr."
end
