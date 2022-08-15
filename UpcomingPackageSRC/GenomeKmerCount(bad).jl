```Has dependencies in old scripts and also is incredibly slow. Just an old proof of concept
```
using BioSequences, FASTX

reader = open(FASTA.Reader, "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna")

function AllReaderKCount(k::Int64, reader::FASTX.FASTA.Reader, kmers::Dict{LongSequence{DNAAlphabet{4}}, Int64})
    #this is the not the fastest version assuming you've already made a kmer dictionary. can be made faster in quite a few ways.
    counts = ReaderKCount(k, reader)
    Counts = copy(counts)
    for key in kmers
        key = first(key)
        if !haskey(counts,key)
            Counts[key] = 0
        end
    end
    return Counts
end


#adding method for k = 1 which has DNA in dictionary. Tbh this probably will never be needed lol
function AllReaderKCount(k::Int64, reader::FASTX.FASTA.Reader, kmers::Dict{DNA, Int64})
    counts = ReaderKCount(k, reader)
    Counts = copy(counts)
    for key in kmers
        key = first(key)
        if !haskey(counts,key)
            Counts[key] = 0
        end
    end
    return Counts
end

#adding a method if you are too lazy to run readerKcount which is actually a pretty bad idea.
function AllReaderKCount(k::Int64, reader::FASTX.FASTA.Reader)
    kmers = genkmersDict(k)
    counts = ReaderKCount(k, reader)
    Counts = copy(counts)
    for key in kmers
        key = first(key)
        if !haskey(counts,key)
            Counts[key] = 0
        end
    end
    return Counts
end
