using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
## available models
include("../src/norm.jl")

p = 1.0
p = .1

data = load_data("30-50", 2016, p, "../data/"; only_positive = true,
                 seed = 1234, opf = true);

out = @time estimate(norm, data, ndc = 4, ngc = 1, normalize = false);

# LP with 15 geo -33208.7, with exp(Gto - Gfrom)
out.out
out.plt[5]
out.plt[6]

