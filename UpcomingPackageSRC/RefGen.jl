#script with functions needed for reference generation. 
"""
genRef(k::Int64,
       reader,
       kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64})

Generate the reference from a FASTA.Reader object of reference sequences and returns a dictionary of the reference.

reader can be a FASTA.Reader object or a string indicating the path of the fasta file.

Is the reference generation algorithm used in the GMA.

can heavily optimized but it doesn't matter too much atm but it will probably not scale too well for very long databases
"""
function genRef(k::Int64, reader::FASTX.FASTA.Reader, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64}) #; returnDict::Bool = true
    len = recordCount(reader)
    answer = Dict{LongSequence{DNAAlphabet{4}}, Float64}()
    for key in kmerDict
        answer[first(key)] = 0.0
    end
    for record in reader
        kmerFreq!(6,FASTA.sequence(record),answer,kmerDict)
    end
    for key in answer
        answer[first(key)] /= len
    end
    return answer
end

function genRef(k::Int64, path::String, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64})
    reader = open(FASTA.Reader,path)
    genRef(k,reader,kmerDict)
end

export genRef

"""
findthr(refseqs::Union{FASTX.FASTA.Reader, String},
        refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
        KD::Dict{LongSequence{DNAAlphabet{4}}, Int64};
        buff::Union{Int64,Float64} = 25)

Prediction of average SED and suggesting a threshold, assuming most indexes do not match.

Its extremely simple and just adds to a the SED of the first reference sequence's KFV to the actual reference KFV
"""
function findthr(refseqs::Union{FASTX.FASTA.Reader, String}, refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}; buff::Union{Int64,Float64} = 25)
    if typeof(refseqs) == String
        refseqs = open(FASTA.Reader, refseqs)
    end
    seq = FASTA.sequence(first(refseqs))
    return Distances.sqeuclidean(kmerFreq(length(first(first(KD))),seq,KD),
    kfv(refKFV,KD)) + buff
end

export findthr

"""
version that scans throug hthe entire reference to see the average SED for the most accurate avg SED but probably doesnt make much of a difference.
"""
function findavgthr(refseqs::Union{FASTX.FASTA.Reader, String}, refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}; buff::Union{Int64,Float64} = 25)
    if typeof(refseqs) == String
        refseqs = open(FASTA.Reader, refseqs)
    end
    SED = []
    KFV = kfv(refKFV,KD)
    k = length(first(first(KD)))
    for record in refseqs
        seq = FASTA.sequence(record)
        push!(SED, Distances.sqeuclidean(kmerFreq(k,seq,KD),KFV))
    end
    return (sum(SED)/recordCount(refseqs)) + buff
end

export findavgthr

AlpacaV = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/AlpacaV.fasta"

lol = genRef(6,AlpacaV,genKmers(6))
#refPlot(lol)

#findthr(AlpacaV,lol,genKmers(6)) #375.797619047619
#findavgthr(AlpacaV,lol,genKmers(6))

#remember to run ApproxMatch.jl and SimpleExplore.jl first.
sixMerDict = genKmers(6)
V3Ref = genRef(6,AlpacaV,sixMerDict)
#V3vec = dictToVec(IMGTRef)
#plot(V3vec)

V3NRef = genRef(6,AlpacaV,genKmers(6,withN=true)) #new genkmer with N
#refPlot(V3NRef)

##V3 reference
V3r = open(FASTA.Reader,V3)

V3ref = genRef(6,V3,sixMerDict)
V3refvec = dictToVec(V3ref)
#refPlot(6,V3refvec)
