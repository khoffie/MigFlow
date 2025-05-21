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
include("../src/choiceset.jl")
include("../src/fullmodel.jl")
include("../src/fullmodel2.jl")
include("../src/othermodels.jl")
include("../src/norm.jl")

data = load_data("30-50", 2017, 0.1, "../data/"; only_positive = true);

out1 = @time estimate(norm, data; norm = "none", type = "conditional");
out1 = @time estimate(norm, data; norm = "both", type = "conditional");
out1 = @time estimate(norm, data; norm = "origin", type = "conditional");
out1 = @time estimate(norm, data; norm = "destination", type = "conditional");

out1 = @time estimate(norm, data; norm = "none", type = "joint");
out1 = @time estimate(norm, data; norm = "both", type = "joint");
out1 = @time estimate(norm, data; norm = "origin", type = "joint");
out1 = @time estimate(norm, data; norm = "destination", type = "joint");

out1 = @time estimate(norm, data; norm = false, type = "joint");
