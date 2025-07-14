#| label: packages-scripts
using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using BenchmarkTools, IterTools, StatProfilerHTML, ReverseDiff
using Suppressor, Distributions
## using Enzyme

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/diagchain.jl")
include("../src/diagrbf.jl")
include("../src/norm.jl")
include("../src/rbf.jl")
include("../src/analysis.jl")

function gendata(ages)
    years = Symbol.("y" .* string.(vcat(2000:2002, 2004:2017)))
    files = readdir("./output"; join = true)

    data = (; (
        age => (; zip(years, reorder(deserialize(file)))...)
        for (age, file) in zip(ages, files)
            )...)
    return data
end
ages = Symbol.("age" .* ["18to25", "25to30", "30to50", "50to65", "above65", "below18"]);

data = gendata(ages);
r = data.agebelow18.y2000
returned(r.mdl.mdl, r.chn)
