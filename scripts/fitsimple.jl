using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, LaTeXStrings, Loess, ADTypes, KernelDensity,
    IterTools, GeoStats, GeoIO, CairoMakie, DynamicPPL, Serialization

include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/modelutils.jl")
include("../src/utils.jl")

include("../models/baseflow.jl")
include("../models/fundamental.jl")
include("../models/gravity.jl")
include("../models/norm.jl")

outp = "./output/"

function defmodel(m, age, year)
    return m(load_data(age, year, 1.0, "../data/"; only_positive = true));
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
years = vcat(2000:2002, 2004:2017)

for  m in [fundamental, norm, gravity]
    for a in ages
        name = "$(m)$(a)"
        results = Vector{EstimationResult}(undef, length(years))
        Threads.@threads for i in eachindex(years)
            results[i] = estimate(defmodel(m, a, years[i]))
        end
        serialize(joinpath(outp, name), results)
        println("$name done")
    end
end
