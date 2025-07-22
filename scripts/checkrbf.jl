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
include("../src/model_helpers.jl")

### Check density RBF
vals = range(-1, 1, 1000)
cx = collect(range(-1, 1, 4))
cy = collect(range(-1, 1, 4))

w = [0, 1, 0, 0,
     0, 0, 0, 0,
     0, 0, 2, 0,
     1, 0, 0, 0]

extract_coefs(out.chn, "ζ")

coefmat(w)
mat, p = plotdensrbf_(Float64.(w), cx, cy, .5, Float64.(vals), "test", nothing);
p

# mdl = norm(load_data("below18", 2014, 1.0, "../data/";
#                      only_positive = true,
#                      seed = 1234, opf = false),
#            normalize = false, ndc = 9, ngcx = 2);
# inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
# @time out = estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
# ## plotdensrbf is in src/diagrbf.jl and basically extracts data from
# ## out and calls plotdensrbf_
# m, p = plotdensrbf(out, nothing);
# m, p = plotgeorbf(out, (-.3, .3));

# p
# out.chn
