```TO DO
    -> make genRefDict version
    -> reader construction.

    IMPORTANT: the plotting of the reference works but,
    something with the function is wrong. the reference generated doesnt have anough counts.

    So, THE REFERENCE GENERATION DOES NOT WORK AT THE MOMENT
```

using BioSequences, FASTX, Plots, Distances

#note: reference generation has dependencies to the other scripts.

```
Pre-generation exploratory analysis
1. N percent
2. average length, Sd, and plot
(Im gonna try combinding these functions into 1 step and have Ns plotted on the bar.)
```
#N percentage. You should inspect the quality of the reference before using it. Ideally as close to 0 as possible.
function percentN(seq::LongSequence{DNAAlphabet{4}})
    allbp = length(seq)
    ncount = 0
    start = 1
    rg = findfirst(ExactSearchQuery(dna"N"), view(seq, start: length(seq)))
    while !isnothing(rg)
        ncount +=1
        start += last(rg)
        rg = findfirst(ExactSearchQuery(dna"N"), view(seq, start: length(seq)))
    end
    return ncount/allbp
end

function percentN(seq::FASTX.FASTA.Record)
    seq = FASTA.sequence(seq)
    percentN(seq)
end

function percentN(collection::FASTX.FASTA.Reader)
    #remember to run the length in ExactMatch first!
    allbp = length(collection)
    ncount = 0
    for seq in collection
        seq = FASTA.sequence(seq)
        start = 1
        rg = findfirst(ExactSearchQuery(dna"N"), view(seq, start: length(seq)))
        while !isnothing(rg)
            ncount +=1
            start += last(rg)
            rg = findfirst(ExactSearchQuery(dna"N"), view(seq, start: length(seq)))
        end
    end
    return ncount/allbp
    close(reader)
end

function percentN(path::String)
    reader = open(FASTA.Reader, path)
    percentN(reader)
end

#In the future I can create fake test data to test this function but its pretty simple and I think it works.
#the following functions are from SimpleExplore.jl

#avfRecLen is from simpleExplore.jl

function genRef(k::Int64, path::String, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64})
    reader = open(FASTA.Reader,path)
    answer = Dict{LongSequence{DNAAlphabet{4}}, Float64}()
    for key in kmerDict
        kmer = first(key)
        answer[kmer] = 0.0
    end
    for record in reader
        rect = kmerFreq(6,FASTA.sequence(record),kmerDict,true)
        for key in answer
            answer[first(key)] += rect[first(key)]
        end
    end
    len = recordCount(open(FASTA.Reader,path))
    for key in answer
        answer[first(key)] /= len
    end
    return answer
end

AlpacaV = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/AlpacaV.fasta"

sixMerDict = genKmers(6)

#remember to run ApproxMatch.jl and SimpleExplore.jl first.
@time IMGTRef = genRef(6,AlpacaV,sixMerDict)
IMGTvec = dictToVec(IMGTRef)

sixMerDict = genKmers(6)
@time V3Ref = genRef(6,AlpacaV,sixMerDict)
V3vec = dictToVec(IMGTRef)
plot(V3vec)

plot(IMGTvec,label = nothing)
title!("Average 6-mer frequency from IMGT Alpaca V genes")
xlabel!("All Unique 6-mers")
ylabel!("Average 6mer count from IMGT reference")
#IT WORKED!!!!!

function refPlot(k::Int64, referece::Vector{Float64},nonzero::Bool = false)
    if nonzero == false
        plot(collect(1:1:length(referece)),referece,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("Average counts")
        title!(string("Average ",k,"-mer distribution in the reference"))
    elseif nonzero == true
        v = Int64[]
        for i in referece
            if i == 0.0
                push!(v,0)
            else
                push!(v,1)
            end
        end
        scatter(v,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("occurence")
        title!(string("all ",k,"-mer occurences in the average of the reference"))

    end
end

refPlot(6,IMGTvec,true)
refPlot(6,IMGTvec)

#THIS DOESNT WORK BUT THE PREVIOUS ONE WORKED.
#I need to make a dict version thats also fast!

#To test I need to construct reader object lol.

#I need to learn how to explort data!!!
close(reader)

##V3 reference

V3r = open(FASTA.Reader,V3)

recordCount(V3)

V3ref = genRef(6,V3,sixMerDict)
V3refvec = dictToVec(V3ref)
avgRecLen(V3r)
refPlot(6, V3refvec)

##

A81 = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/genBank/A81.fasta"
A91 = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/genBank/A91.fasta"
J71 = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/genBank/J71.fasta"
J81 = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/genBank/J81.fasta"
