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
include("../src/chebies.jl")
include("../src/diagchain.jl")
include("../src/diagrbf.jl")
include("../src/norm.jl")
include("../src/rbf.jl")

corround(x, y) = round(cor(x, y), digits = 2)

p = 1.0
p = 0.1
data = load_data("18-25", 2017, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);

AD = AutoForwardDiff()
## AD = AutoReverseDiff()
## AD = AutoEnzyme()
mdl = norm(data;
           kdens = 2.0,
           kgeo = 2.0,
           ndc = 4, ngc = 4,
           normalize = false);

Random.seed!(1234)
chn = Turing.sample(mdl.mdl, NUTS(100, .7; adtype = AD), 1, progress = true)
out = chaindiag(chn, mdl, data);

# plot(out.plts[1:4]...)
# out.plts[5]
# out.plts[6]
# chn

## tree depth 4 or 5
## check how NUTS works
## rescale y axis, such that 1m = 1m

# out = @time estimate(norm, data;
#                      model_kwargs = (; ndc = 4, ngc = 4, normalize = false),
#                      optim_kwargs = (; ad = AD));
