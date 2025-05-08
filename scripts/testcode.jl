using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings
using BenchmarkTools, ADTypes, ReverseDiff, PrettyTables

include("../src/utils.jl")
include("../src/othermodels.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/fullmodel.jl")
include("../src/chebies.jl")
include("../src/choiceset.jl")

mdat = gen_mdat(load_data("30-50", 2017, 0.1, "../data/");
                distscale = 100.0, ndc = 28, ngc = 28)

out = @time estimate(choice, mdat; normalize = false);
out2 = @time estimate(choice, mdat; normalize = true);
out.out
out2.out
