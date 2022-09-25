```TO DO
    -> make genRefDict version
    -> reader construction.

    IMPORTANT: the plotting of the reference works but,
    something with the function is wrong. the reference generated doesnt have anough counts.

    So, THE REFERENCE GENERATION DOES NOT WORK AT THE MOMENT
```

using BioSequences, FASTX, Plots, Distances

#note: reference generation has dependencies to the other sripts.
#In the future I can create fake test data to test this function but its pretty simple and I think it works.
#the following functions are from SimpleExplore.jl

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

@time V3NRef = genRef(6,AlpacaV,genKmers(6,withN=true)) #new genkmer with N
refPlot(V3NRef)

plot(IMGTvec,label = nothing)
title!("Average 6-mer frequency from IMGT Alpaca V genes")
xlabel!("All Unique 6-mers")
ylabel!("Average 6mer count from IMGT reference")
#IT WORKED!!!!!

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
