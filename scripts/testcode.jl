using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings
using ADTypes, ReverseDiff, Loess

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
include("../src/gen_mdat.jl")
## available models
include("../src/norm.jl")

p = 1.0
p = .1

data = load_data("30-50", 2017, p,
                 "../data/";
                 only_positive = true,
                 seed = 1234,
                 opf = false);

out1 = @time estimate(norm, data; ndc = 1, ngc = 1, normalize = true);
out2 = @time estimate(norm, data; ndc = 1, ngc = 1, normalize = false);

out1.plt[5]
