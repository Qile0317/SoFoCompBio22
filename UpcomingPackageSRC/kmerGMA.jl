module kmerGMA
    using Plots, BioSequences, FASTX, WebBlast, Distances, Statistics,
    ProgressMeter, StringDistances, StatsBase, StatsPlots, DataFrames
    include("simpleExplore.jl")
    include("ExactMatch.jl")
    include("RefGen.jl")
    include("GMA.jl")
    #include("GMA_Nless.jl")
end
