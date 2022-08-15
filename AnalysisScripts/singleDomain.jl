``` just a super messy script used to find simgle domain nanobodies from amino acid
sequences based on the 4 hallmarks at AA positions 42,49,50,52 from this paper.
```

using BioSequences,FASTX

Dpath = "C:/Users/lu_41/Desktop/Sofo Prok/AA.fasta"

reader = open(FASTA.Reader,Dpath)

#nanobody
SDAB = Dict{String,Int64}()
for record in reader
    seq = FASTA.sequence(record)
    score = 0
    if seq[42] == AA_Y || seq[42] == AA_F
        score += 1
    end
    if seq[49] == AA_E || seq[49] == AA_V
        score+=1
    end
    if seq[50] == AA_G || seq[50] == AA_R
        score+=1
    end
    if seq[52] == AA_G
        score += 1
    end
    SDAB[FASTA.identifier(record)] = score
end
SDAB

#changing string name
nd = Dict{String,Int64}()
for key in SDAB
    iden = first(key)
    fam = iden[10:end-8]
    final ="Vicpac "*fam
    nd[final] = SDAB[iden]
end
return nd

refe = ["Vicpac IGHV1-1*01 F", "Vicpac IGHV1S3*01 F", "Vicpac IGHV1S4*01 F", "Vicpac IGHV1S5*01 P", "Vicpac IGHV1S6*01 F", "Vicpac IGHV3-1*01 F", "Vicpac IGHV3-2*01 F", "Vicpac IGHV3-3*01 F", "Vicpac IGHV3S1*01 F", "Vicpac IGHV3S10*01 F", "Vicpac IGHV3S11*01 F", "Vicpac IGHV3S12*01 F", "Vicpac IGHV3S15*01 F", "Vicpac IGHV3S16*01 F", "Vicpac IGHV3S19*01 P", "Vicpac IGHV3S2*01 F", "Vicpac IGHV3S25*01 F", "Vicpac IGHV3S26*01 F", "Vicpac IGHV3S28*01 F", "Vicpac IGHV3S29*01 F", "Vicpac IGHV3S30*01 F", "Vicpac IGHV3S31*01 F", "Vicpac IGHV3S32*01 F", "Vicpac IGHV3S33*01 F", "Vicpac IGHV3S34*01 F", "Vicpac IGHV3S35*01 F", "Vicpac IGHV3S36*01 F", "Vicpac IGHV3S37*01 F", "Vicpac IGHV3S39*01 F", "Vicpac IGHV3S40*01 F", "Vicpac IGHV3S41*01 F", "Vicpac IGHV3S42*01 F", "Vicpac IGHV3S44*01 F", "Vicpac IGHV3S46*01 F", "Vicpac IGHV3S53*01 F", "Vicpac IGHV3S54*01 F", "Vicpac IGHV3S55*01 F", "Vicpac IGHV3S57*01 F", "Vicpac IGHV3S58*01 F", "Vicpac IGHV3S59*01 F", "Vicpac IGHV3S6*01 F", "Vicpac IGHV3S60*01 F", "Vicpac IGHV3S61*01 F", "Vicpac IGHV3S62*01 F", "Vicpac IGHV3S63*01 F", "Vicpac IGHV3S64*01 F", "Vicpac IGHV3S65*01 F", "Vicpac IGHV3S66*01 F", "Vicpac IGHV3S67*01 F", "Vicpac IGHV3S68*01 F", "Vicpac IGHV3S9*01 F", "Vicpac IGHV4S1*01 F", "Vicpac IGHV4S11*01 P", "Vicpac IGHV4S2*01 F", "Vicpac IGHV4S3*01 F", "Vicpac IGHV4S5*01 F", "Vicpac IGHV4S6*01 F", "Vicpac IGHV4S7*01 F", "Vicpac IGHV4S8*01 F", "Vicpac IGHV4S9*01 P"]

nref = []
for i in refe
    push!(nref,i[1:end-2])
end
nref

for i in nref
    prev = nothing
    if i!=prev
        prev = i
    elseif i == prev
        return true
    end
end

answe = []
for i in nref
    if haskey(nd,i)
        push!(answe,nd[i])
    end
end

answe

println(answe)

ANSWE = Dict()
for i in nref
    for key in nd
        ident = first(key)
        if i == ident
            ANSWE[ident] = nd[ident]
        end
    end
end

ANSWE
