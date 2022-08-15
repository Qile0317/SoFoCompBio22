#script with me loading in data from IMGT analyss of seqeuenced data

using CSV, DataFrames, Statistics, BioSequences, Plots, ProgressMeter, StringDistances, SeqUMAP, StatsBase, DelimitedFiles, FASTX
#This happened when installing NextGenSeqUtils
```
✗ NextGenSeqUtils
10 dependencies successfully precompiled in 84 seconds (231 already precompiled)
1 dependency errored. To see a full report either run `import Pkg; Pkg.precompile()` or load the package
```

#IMPORTANT::

```
each isotype file will have the same primer sequence (with some error tolerance) because thats how I split the sequencing file into isotypes
2:28
and you know the primer we used for each isotype
```

#loading data

IgE = "C:/Users/lu_41/Desktop/Sofo Prok/__MACOSX/IgE.tsv"
IgG1 = "C:/Users/lu_41/Desktop/Sofo Prok/__MACOSX/IgG1.tsv"
IgG2b = "C:/Users/lu_41/Desktop/Sofo Prok/__MACOSX/IgG2b.tsv"
IgG2c = "C:/Users/lu_41/Desktop/Sofo Prok/__MACOSX/IgG2c.tsv"
IgM_A = "C:/Users/lu_41/Desktop/Sofo Prok/__MACOSX/IgM_A.tsv"
IgM_B = "C:/Users/lu_41/Desktop/Sofo Prok/__MACOSX/IgM_B.tsv"
IgA = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/IMGT_IgA/aair/vquest_airr.tsv"

IgEdf = CSV.read(IgE,DataFrame,delim = '\t')
IgG1df = CSV.read(IgG1,DataFrame,delim = '\t')
IgG2bdf = CSV.read(IgG2b,DataFrame,delim = '\t')
IgG2cdf = CSV.read(IgG2c,DataFrame,delim = '\t')
IgM_Adf = CSV.read(IgM_A,DataFrame,delim = '\t')
IgM_Bdf = CSV.read(IgM_B,DataFrame,delim = '\t')
IgAdf = CSV.read(IgA,DataFrame,delim = '\t')

function oppo(dat::Vector)
    ansvec= Float64[]
    for i in dat
        if !Base.ismissing(i)
            push!(ansvec,100-i)
        end
    end
    return ansvec
end

Eshm = oppo(IgEdf[:,:v_identity])
G1shm = oppo(IgG1df[:,:v_identity])
MAshm = oppo(IgM_Adf[:,:v_identity])
MBshm = oppo(IgM_Bdf[:,:v_identity])
G2bshm = oppo(IgG2bdf[:,:v_identity])
G2cshm = oppo(IgG2cdf[:,:v_identity])
Ashm = oppo(IgAdf[:,:v_identity])
```
#conversion to biosequence
function c2dna(vec::Vector)
    nseqs = []
    for i in vec
        push!(nseqs,LongSequence{DNAAlphabet{4}}(i))
    end
    return nseqs
end
```
#getting the seqeucnes
Es = IgEdf[:, :sequence]
G1s = IgG1df[:, :sequence]
G2bs = IgG2bdf[:, :sequence]
G2cs = IgG2cdf[:, :sequence]
MAs = IgM_Adf[:, :sequence]
MBs = IgM_Bdf[:, :sequence]
As = IgAdf[:, :sequence] #i also have the original files

Alls = [Es,G1s,G2bs,G2cs,MBs,As]
AllS = [Es,G1s,G2bs,G2cs,MAs,MBs,As]
Allnames = ["IgE","IgG1", "IgG2b","IgG2c","IgM_A", "IgM_B","IgA"]

Ev = "CCCAGAAACCAACCGTCTTCCCCTTGACCTGTTGCAAAAACACCACCG"

G1v = "AGGCCCCATCGGTCTATCCTCTGACTGCTAGATGCGGGGACACGC"

G2bv = "ACCACAACCACAACCACAACCCAATCCTACAACAGAATCCAAGTGTCCCAAATGTCC"

G2cv = "GCACCACAGCGAAGACCCCAGCTCCAAGTGTCCCAAATGCCCAGG"

MBv = "TGGCCTTGGGCTGCCTAGCCCGGGACTTCCTGCCTGGCTCCATC"

Av = "CCAGCCCCAGCGTCTTCCCGCTGGGCCCCAGCTATGACAAGGCATC"

Allv = [Ev,G1v,G2bv,G2cv,MBv,Av]

function checkseq(seqs::Vector{String},query::String)
    query = lowercase(query)
    output = Int64[]
    l = length(query)
    for seq in seqs
        push!(output,Levenshtein()(query,seq[end-l+1:end]))
    end
    return output
end

#making the violin plots for the identifier region to the right region
Plots.violin(["IgE"],checkseq(Es,Ev),label=nothing)
Plots.violin!(["IgG1"],checkseq(G1s,G1v),label=nothing)
Plots.violin!(["IgG2b"],checkseq(G2bs,G2bv),label=nothing)
Plots.violin!(["IgG2c"],checkseq(G2cs,G2cv),label=nothing)
Plots.violin!(["IgM_B"],checkseq(MBs,MBv),label=nothing)
Plots.violin!(["IgA"],checkseq(As,Av),label=nothing)
xlabel!("Filtered Isotypes")
ylabel!("Levenshtein distance to identifier region")

function porpBelow(ldv::Vector{Int64},thr::Int64 = 3)
    t = 0
    f = length(ldv)
    for i in ldv
        if i < thr
            t+= 1
        end
    end
    final = string(t*100/f)
    println(final[1:5]*"%")
end

porpBelow(checkseq(As,Av))

function allcheckseq(seqs::Vector{Vector{String}},query::Vector{String})
    query = [lowercase(i) for i in query]
    output = Int64[]
    l = [length(i) for i in query]
    for i in 1:length(seqs)
        seqVector = seqs[i]
        for seq in seqVector
            push!(output,Levenshtein()(query[i],seq[end-l[i]+1:end]))
        end
    end
    return output
end

AllDist = allcheckseq(Alls,Allv)
histogram(AllDist, label = nothing, bins = 30)
xlabel!("Levenshtein distance fo identifier region")
ylabel!("amount of sequences")
title!("Histogram")



#making the violin plot of lengths

function aLen(seqVec::Vector{Vector{String}})
    outvec = Vector{Int64}[]
    for seqs in seqVec
        lenVec = Int64[]
        for seq in seqs
            push!(lenVec,length(seq))
        end
        push!(outvec,lenVec)
    end
    return outvec
end

AllLen = aLen(AllS)

#violin plot. i tried making a function for this but it wouldnt make any plot but still runs.
Plots.violin(["IgE"],AllLen[1],label=nothing)
Plots.violin!(["IgG1"],AllLen[2],label=nothing)
Plots.violin!(["IgG2b"],AllLen[3],label=nothing)
Plots.violin!(["IgG2c"],AllLen[4],label=nothing)
Plots.violin!(["IgM_A"],AllLen[5],label=nothing)
Plots.violin!(["IgM_B"],AllLen[6],label=nothing)
Plots.violin!(["IgA"],AllLen[7],label=nothing)
xlabel!("Filtered Isotypes")
ylabel!("Merged read length")

#verifying violin plot by looking at primer to J region:
```
IgE: 26
IgG1: 30
IgG2b: 47
IgG2c: 22
IgM_B: 90
IgA: 30

```

#getting the alignments
E_Al = IgEdf[:, :sequence_alignment]
G1_Al = IgG1df[:, :sequence_alignment]
G2b_Al = IgG2bdf[:, :sequence_alignment]
G2c_Al = IgG2cdf[:, :sequence_alignment]
MA_Al = IgM_Adf[:, :sequence_alignment]
MB_Al = IgM_Bdf[:, :sequence_alignment]
A_Al = IgAdf[:, :sequence_alignment]

#removing dots from alignment, taking first bp basepairs and making a countmap
function unAlign(Alignment::Vector, bp::Int64, cm::Bool)
    answer = []
    @showprogress for seq in Alignment
        if ismissing(seq) == false
            currseq = ""
            currlen = 0
            for nt in seq
                if currlen < bp
                    if nt != '.'
                        currseq = currseq*nt
                        currlen += 1
                    end
                end
            end
            push!(answer,currseq)
        else
            push!(answer,missing)
        end
    end
    if cm
        return countmap(answer)
    else
        return answer
    end
end

function unAlign(Alignment::Vector, bp::Int64)
    answer = String[]
    @showprogress for seq in Alignment
        if ismissing(seq) == false
            currseq = ""
            currlen = 0
            for nt in seq
                if currlen < bp
                    if nt != '.'
                        currseq = currseq*nt
                        currlen += 1
                    end
                end
            end
            push!(answer,currseq)
        end
    end
    return answer
end

#unaligning and not countmapping
Ensa = unAlign(E_Al, 285)
G1nsa = unAlign(G1_Al, 285)
G2bnsa = unAlign(G2b_Al,285)
G2cnsa = unAlign(G2c_Al,285)
MAnsa = unAlign(MA_Al,285)
MBnsa = unAlign(MB_Al,285)
Ansa = unAlign(A_Al,285)

#unaligning and countmapping (note, cm is wrong here)
Ens = unAlign(E_Al, 285)
G1ns = unAlign(G1_Al, 285)
G2bns = unAlign(G2b_Al,285)
G2cns = unAlign(G2c_Al,285)
@time MAns = unAlign(MA_Al,285) # 15.579581 seconds (297.33 M allocations: 21.294 GiB, 45.18% gc time)
MBns = unAlign(MB_Al,285)
Ans = unAlign(A_Al,285)

Alls = [Ens,G1ns,G2bns,G2cns,MAns,MBns,Ans]
AllsNames = ["IgE","IgG1x","IgG2b","IgG2c","IgM_A","IgM_B","IgA"]

#finding the top 100. surprisingly fast. Isnt really needed as im just printing to REPL
function topCounts(ct::Dict{Any,Int64}, thr::Int64, countDict::Bool = false)
    sct =  reverse(sort(collect(zip(values(ct),keys(ct)))))
    if countDict == false
        answer = String[]
        curr = 0
        for t in sct
            if curr < thr
                if !ismissing(last(t))
                    push!(answer, last(t))
                    curr += 1
                end
            end
        end
        return answer
    else
        answer = Dict{String,Int64}()
        curr = 0
        for t in sct
            if curr < thr
                if !ismissing(last(t))
                    curr += 1
                    answer[string(name)*" seq_"*string(curr)] = first(t)
                end
            end
        end
        return answer
    end
end

thE = topCounts(Ens,100)
thG1 = topCounts(G1ns,100)
thG2b = topCounts(G2bns,100)
thG2c = topCounts(G2cns,100)
thMA = topCounts(MAns,100)
thMB = topCounts(MBns,100)
thA = topCounts(Ans,100)

#in the futurue ill get more familiar with creating and writing into files but for now i can just print
#but this version doesnt do all of them cumulatively. So im gonna make another ver
function printTC(ct::Dict{Any,Int64}, thr::Int64, name::String)
    sct =  reverse(sort(collect(zip(values(ct),keys(ct)))))
    curr = 0
    for t in sct
        if curr < thr
            if !ismissing(last(t))
                curr += 1
                println(">"*name*" seq "*string(curr)*" | count = "*string(first(t)))
                println(last(t))
            end
        end
    end
end


printTC(Ans,100,"IgA")

ct = MAns
nct = reverse(sort(collect(zip(values(ct),keys(ct))))) #1435, idk why it doesnt show up in repl
println(">IgM_A seq 1 | count = 1435") #done. the seqs are in a fasta file now

function printTC(ct::Vector{Dict{Any,Int64}}, thr::Int64, names::Vector{String})
    globalcurr = 0
    dcount = 0
    for dict in ct
        curr = 0
        dcount += 1
        sct =  reverse(sort(collect(zip(values(dict),keys(dict)))))
        for t in sct
            if curr < thr
                if !ismissing(last(t))
                    curr += 1
                    globalcurr += 1
                    println(">"*names[dcount]*" seq "*string(curr)*" | seq_"*string(globalcurr)*" | count = "*string(first(t)))
                    println(last(t))
                end
            end
        end
    end
end

printTC(Alls,100,AllsNames)

#all previous functions were for the phylogeny

#doing a sequmap to check similarity and it looks cool
#if this doesnt work ill come back to it. (hasnt worked yet, also havent figured out labelling)

UAs = String[]
for i in As
    push!(UAs, uppercase(i))
end

@time proj = sequmap(UAs, 2; k = 6, min_dist = 1e-3)
#idk how to plot the matrix. obviously regular plots didnt work.

#I gotta try make a circular phylogenetic tree as well in fastree and figtree. check the linux pipeline


```
you probably want to look at sequence alignments and pick out a little signature sequence for each isotype, right after the primer
and then match that to the sequence, right after the primer
allowing for some allelic variation, so it mustn't just be perfect matching
```
```
#the enxt part uses a function from ExactMatch.jl which is findall. Here i added a method without needing an empty vector. Problem I need percentage mismatches.
function ApproxFA(q::ApproximateSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, seq::LongSequence{DNAAlphabet{4}}, mismatch::Int64)
    answer= UnitRange[]
    start = 1
    rg = findfirst(q, mismatch, view(seq, start: length(seq)))
    while !isnothing(rg)
        push!(answer, start-1+first(rg): start-1+last(rg))
        start += last(rg)
        rg = findfirst(q, mismatch, view(seq, start: length(seq)))
    end
    if answer == UnitRange[]
        return nothing
    else
        return answer
    end
end

ApproxFA(E_siga,LongSequence{DNAAlphabet{4}}(Ep["seq_1"]), 3)

Ep["seq_1"][502:525]

function d2s(seq::LongSequence{DNAAlphabet{4}})
    for i in seq
        seq

function ApproxFF(q::String, seqs::Dict{}, mismatch::Int64)
    lq = length(q)
    qu = ApproximateSearchQuery(LongSequence{DNAAlphabet{4}}(q))
    resultD = Dict{String, Any}() #unitrange?
    for key in seqs
        loc = findfirst(qu,mismatch,LongSequence{DNAAlphabet{4}}(last(key)))
        if loc != nothing
            m = last(key)
            resultD[first(key)] = [loc,m[loc],Levenshtein()(m[loc],q)]
        elseif loc == nothing
            resultD[first(key)] = "no match"
        end
    end
    return resultD
end

E = ApproxFF(E_sig, Ep, 2)

G1 = ApproxFF(G1_sig,G1p,2)
```
##Screw the prev function. I did this for all of them except IgM A to check them in aliview

function printall(Ep::Dict, name::String, am::Int64 = 200)
    co = 0
    for key in Ep
        if co < am+1
            println(">"*name)
            println(last(key))
            println("")
            co += 1
        end
    end
end

function printall(Ep::Dict, name::String)
    for key in Ep
        println(">"*name)
        println(last(key))
        println("")
    end
end

printall(G2bp, "IgG2B")

##
```
function c2ds(seq::String)
    seq = ApproximateSearchQuery(LongSequence{DNAAlphabet{4}}(seq))
    return seq
end
```
#EXTREMELY IMPORTANT
Ev = "CCCAGAAACCAACCGTCTTCCCCTTGACCTGTTGCAAAAACACCACCG"

G1v = "AGGCCCCATCGGTCTATCCTCTGACTGCTAGATGCGGGGACACGC"

G2bv = "ACCACAACCACAACCACAACCCAATCCTACAACAGAATCCAAGTGTCCCAAATGTCC"

G2cv = "GCACCACAGCGAAGACCCCAGCTCCAAGTGTCCCAAATGCCCAGG"

MBv = "TGGCCTTGGGCTGCCTAGCCCGGGACTTCCTGCCTGGCTCCATC"

Av = "CCAGCCCCAGCGTCTTCCCGCTGGGCCCCAGCTATGACAAGGCATC"

function checkseq(seqs::Vector{String},query::String)
    query = lowercase(query)
    output = Int64[]
    l = length(query)
    for key in seqs
        push!(output,Levenshtein()(query,last(key)[end-l+1:end]))
    end
    return output
end


function csz(seqs::Dict{String,String},query::String)
    query = lowercase(query)
    output = 0
    l = length(query)
    for key in seqs
        if Levenshtein()(query,last(key)[end-l+1:end]) == 0
            output += 1
        end
    end
    return output
end

csz(Ep,Ev) + csz(G1p,G1v) + csz(G2bp,G2bv) + csz(G2cp,G2cv) + csz(MBp,MBv) + csz(Ap,Av)

function porpNZ(d::Dict{String,Int64})
    z = 0
    for key in d
        if last(key) == 0
            z += 1
        end
    end
    return 1-(z/length(d))
end

function plotMat(d::Dict{String,Int64},title::String)
    vec = Int64[]
    for key in d
        push!(vec,last(key))
    end
    Plots.histogram(vec,label=nothing, titlefontsize= 10, linewidth=0.1)
    title!(title*" mismatch distribution")
    xlabel!("levenshtein distance to identifier region")
    ylabel!("number of sequences")
end

Eplt = plotMat(Emat,"IgE")
G1plt = plotMat(G1mat,"IgG1")
G2bplt = plotMat(G2bmat,"IgG2b")
G2cplt = plotMat(G2cmat,"IgG2c")
MBplt = plotMat(MBmat,"IgM_B")
#not IgM_A remember!!
Aplt = plotMat(Amat, "IgA")

plot(Eplt,G1plt,G2bplt,G2cplt,MBplt, Aplt, legend = false) #doenst work lol

function notZ(d::Dict{String,Int64},threshold::Int64) #gives the mismatched sequences
    nd = Dict{String,Int64}()
    for key in d
        if last(key) >= threshold
            nd[first(key)] = d[first(key)]
        end
    end
    return nd
end

#here are the sequences that were not perfect matches.
Enz = notZ(Emat,10)
G1nz = notZ(G1mat,10)
G2bnz = notZ(G2bmat,10)
G2cnz = notZ(G2cmat,10)
MBnz = notZ(MBmat,10)
Anz = notZ(Amat,10)

function id2s(dist::Dict{String, Int64},re::Dict{String, String},name::String)
    a = Dict{String,String}()
    for key in dist
        ide = first(key)
        seq = re[ide]
        a[ide] = seq
    end
    printall(a,name)
end

id2s(Enz,Ep,"IgE mismatch | thr=10")

id2s(G1nz,G1p,"IgG1x mismatch | thr=10")

id2s(G2bnz,G2bp,"IgG2b mismatch | thr=10")

id2s(G2cnz,G2cp,"IgG2c mismatch | thr=10")

id2s(MBnz,MBp, "IgM_B mismatch | thr=10")

id2s(Anz,Ap, "IgA mismatch | thr=10")
