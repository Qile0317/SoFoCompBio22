```
Problem: Idk what happened but I think GenomeMatch is brokeen :/
```

##import packages
using BioSequences, FASTX, Plots
import Base.findall, Base.length

##
reader = open(FASTA.Reader, "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna")
path = "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna"
firstr = first(reader)
f = FASTA.sequence(firstr)
q = dna"agtatcactaattatcagagaaatgcaaatcaaaactac"
query = ExactSearchQuery(q)

bseq = dna"CCCCCC"^10

## I can also make it return a dictionary.
function FirstMatch(reader::FASTX.FASTA.Reader{}, query::LongSequence{DNAAlphabet{4}})
    query = ExactSearchQuery(query)
    for record in reader
        find = findfirst(query, FASTA.sequence(LongSequence{DNAAlphabet{4}}, record))
        if !isnothing(find)
            println(find," ", FASTA.identifier(record))
        end
    end
end

@time FirstMatch(reader,q)

##
function FindAll(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, seq::LongSequence{DNAAlphabet{4}}, answer::Vector{UnitRange})
    start = 1
    rg = findfirst(q, view(seq, start: length(seq)))
    while !isnothing(rg)
        push!(answer, start-1+first(rg): start-1+last(rg))
        start += last(rg)
        rg = findfirst(q, view(seq, start: length(seq)))
    end
    if answer == UnitRange[]
        return nothing
    else
        return answer
    end
end

answer = UnitRange[]
@time FindAll(query, f, answer) #0.000727 seconds (3 allocations: 160 bytes)

#overlap
function FindAllOverlap(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, seq::LongSequence{DNAAlphabet{4}}, answer::Vector{UnitRange})
    start = 1
    rg = findfirst(q, view(seq, start: length(seq)))
    while !isnothing(rg)
        push!(answer, start-1+first(rg): start-1+last(rg))
        start += first(rg)
        rg = findfirst(q, view(seq, start: length(seq)))
    end
    if answer == UnitRange[]
        return nothing
    else
        return answer
    end
end

answer = UnitRange[]
@time FindAllOverlap(query,f,answer)

answer = UnitRange[]
@time FindAllOverlap(query,bseq,answer)

function RecordMatch(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, seq::LongSequence{DNAAlphabet{4}}, overlap::Bool)
    answer = UnitRange[]
    if overlap==true
        FindAllOverlap(q, seq, answer)
    elseif overlap == false
        FindAll(q,seq, answer)
    end
end

@time RecordMatch(query,f,true)
@time RecordMatch(query,bseq,true)

##
function GenomeMatch(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, Reader::FASTX.FASTA.Reader{},overlap::Bool = true)
    identify = Dict{String,Vector{UnitRange{Int64}}}()
    for record in reader
        seq = FASTA.sequence(record)
        RM = RecordMatch(q,seq,overlap)
        if !isnothing(RM)
            identify[FASTA.identifier(record)] = RM
        end
    end
    if identify == Dict{String, Vector{UnitRange{Int64}}}
        return "no match"
    else
        return identify
    end
    close(reader)
end

function GenomeMatch(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, path::String, overlap::Bool=true)
    reader = open(FASTA.Reader,path)
    GenomeMatch(q,reader,overlap)
end

@time matchess = GenomeMatch(query, path)

##
#rember:
FASTX.FASTA.seqlen(firstr) # 0.000001 seconds
#is much faster.

#and scaling to the entire assembly (once again by adding to the base method length)
function length(reader::FASTX.FASTA.Reader)
    len = 0
    for record in reader
        len += FASTX>FASTA.seqlen(record)
    end
    return len
end

@time length(reader) #2.747

## function to record and store cumulative lengths of the BEGINNING of each record in a dictionary
function cflength(reader::FASTX.FASTA.Reader{})
    lengthmap = Dict{String,Int64}()
    clength = 0
    prev = 0
    for record in reader
        identifier = FASTA.identifier(record)
        clength += prev
        prev = FASTX.FASTA.seqlen(record)
        lengthmap[identifier] = clength
    end
    return lengthmap
end

clengths = cflength(reader)
clengths[FASTA.identifier(firstr)] #works as intended :)

#plotting - ver needing clength and length(reader) for speed
```
i was thinking to maybe change the data strucutre for the vector of unitranges to something like a hashmap
to have constant lookup times instead of iterating through the vector and making the time complexity O(n^2)
```
function PlotQueryMatches(matches::Dict{String, Vector{UnitRange{Int64}}},  clengths::Dict{String,Int64}, rlength::Int64)
    x = Int64[]
    for match in matches
        currlen = clengths[first(match)]
        c = matches[first(match)]
        for ur in c
            f = first(ur) + currlen
            push!(x,f)
        end
    end
    push!(x,rlength) #for a slightly more elegant graph i could maybe try find a way to extend the x axis to rlength
    y = fill(1,length(x))
    scatter(x,y, title = "query match locations", label = "first position of match")
    xlabel!("matches along genome")
end

@time PlotQueryMatches(matchess,clengths,rlength)

#plotting - slower ver not needing any prereq except matches and reader
function SlowerPQM(matches::Dict{String, Vector{UnitRange{Int64}}},reader::FASTX.FASTA.Reader{})
    clengths = cflength(reader)
    rlength = length(reader)
    PlotQueryMatches(matches,clengths,rlength)
end
#the plotting has to be displayed somehow

@time SlowerPQM(matchess, reader)

## function that does everything by finding matches AND plotting.
function MatchPlot(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, Reader::FASTX.FASTA.Reader{},overlap::Bool)
    matches = GenomeMatch(q,reader,overlap)
    return matches
    return SlowerPQM(matches,reader)
end

MatchPlot(query, reader, true)
