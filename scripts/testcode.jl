using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using BenchmarkTools, IterTools, StatProfilerHTML, ReverseDiff
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

mdl = norm(load_data("18-25", 2014, 1.0, "../data/";
                     only_positive = true,
                     seed = 1234, opf = false),
           normalize = false, ndc = 9, ngcx = 2);
inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
@time out = estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
m, p = plotdensrbf(out, nothing);
p
