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

p = 1.0
p = 0.1
data = load_data("18-25", 2015, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);

AD = AutoForwardDiff()
## AD = AutoReverseDiff()
## AD = AutoEnzyme()
mdl = norm(data,
           kdens = 2.0,
           kgeo = 2.0,
           ndc = 4, ngcx = 2,
           normalize = false);
inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);

out = @time estimate(mdl, optim_kwargs = (; show_trace = false,
                                          inits = inits));
plot(out.plts[1:4]...)
plot(out.plts[5])
plot(out.plts[6])
