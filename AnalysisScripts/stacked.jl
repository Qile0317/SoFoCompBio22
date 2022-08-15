using CSV, DataFrames, Statistics, BioSequences, PyPlot, Plots, StatsPlots, StringDistances, SeqUMAP, StatsBase, DelimitedFiles

iso = "C:/Users/lu_41/Desktop/Sofo Prok/IsoMod.tsv"
SHM = "C:/Users/lu_41/Desktop/Sofo Prok/Iso.tsv"

isodf = CSV.read(SHM,DataFrame)
println(names(isodf))

#ignore
function di(n::Float64)
        return n/153.26
end

function dia(v::Vector)
        nv = Float64[]
        for i in v
                push!(nv,di(i))
        end
        return nv
end

MApe = isodf[:, :IgM_A]
MBpe = isodf[:, :IgM_B]
G1pe = isodf[:, :IgG1]
G2bpe = isodf[:, :IgG2b]
G2cpe = isodf[:, :IgG2c]
Epe = isodf[:, :IgE]

#ignore this part. how do I do this for A?
Apath = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/IMGT_IgA/8_V-REGION-nt-mutation-statistics.txt"
Asum = CSV.read(Apath,DataFrame, delim = '\t')
println(names(Asum))

# In PyPlot backend, if we use chars like 'A':'L', ticks are displayed with "PyWrap".
ticklabel = isodf[:, :Hallmarks]
ticklabel[end] = 0

groupedbar([MApe MBpe G1pe G2bpe G2cpe Epe],
        bar_position = :stack,
        #bar_width=0.7,
        xticks=(1:61, ticklabel),
        label=["IgM_A" "IgM_B" "IgG1" "IgG2b" "IgG2c" "IgE"])

xlabel!("SdAb Hallmark counts from IMGT germline references")
ylabel!("Cumulative isotype counts (%)")

##using data from validtion.jl:
Plots.violin(["IgM_A"],MAshm, label=nothing)
Plots.violin!(["IgM_B"],MBshm, label=nothing)
Plots.violin!(["IgG1"],G1shm, label=nothing)
Plots.violin!(["IgG2b"],G2bshm, label=nothing)
Plots.violin!(["IgG2c"],G2cshm, label=nothing)
Plots.violin!(["IgE"],Eshm, label=nothing)
Plots.violin!(["IgA"],Ashm, label=nothing)
xlabel!("Isotypes")
ylabel!("SHM (%)")


```
#ignore, these are averages
MAshm = isodf[:, :IgM_A_VSHM]
MBshm = isodf[:, :IgM_B_VSHM]
G1shm = isodf[:, :IgG1_VSHM]
G2bshm = isodf[:, :IgG2b_VSHM]
G2cshm = isodf[:, :IgG2c_VSHM]
Eshm = isodf[:, :IgE_VSHM]

Plots.violin(["IgM_A"],MAshm, label=nothing)
Plots.violin!(["IgM_B"],MBshm, label=nothing)
Plots.violin!(["IgG1"],G1shm, label=nothing)
Plots.violin!(["IgG2b"],G2bshm, label=nothing)
Plots.violin!(["IgG2c"],G2cshm, label=nothing)
Plots.violin!(["IgE"],Eshm, label=nothing)
xlabel!("Isotypes")
ylabel!("SHM (%)")

#ignore, wrong plot
groupedbar([MAshm MBshm G1shm G2bshm G2cshm Eshm],
        bar_position = :stack,
        #bar_width=0.7,
        xticks=(1:61, ticklabel),
        label=["IgM_A" "IgM_B" "IgG1" "IgG2b" "IgG2c" "IgE"])

xlabel!("SdAb Hallmark counts from IMGT germline references")
ylabel!("Cumulative SHM (%)")
```
