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
include("../src/norm.jl")
include("../src/fullmodel.jl")
include("../src/fullmodel2.jl")
include("../src/othermodels.jl")

data = load_data("30-50", 2017, 0.1, "../data/"; only_positive = true);
mdat = gen_mdat(data; type = "conditional", distscale = 100.0, ndc = 28, ngc = 1);
out1 = @time estimate(distonly, mdat);
out1.out
out1 = @time estimate(norm, mdat);
out1.out
