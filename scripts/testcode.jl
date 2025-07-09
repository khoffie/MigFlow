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
mdl = norm(data,
           kdens = 2.0,
           kgeo = 2.0,
           ndc = 4, ngcx = 2,
           normalize = false);

out = @time estimate(mdl, optim_kwargs = (; show_trace = false));

rbfinits(N, σ, t = 90.0) = clamp.(rand(MvNormal(zeros(N), σ^2 *I(N))), -t, t)


function initialize(a, ndc, ngcx, ngcy)
    density = rbfinits(ndc, 40.0)
    geo = rbfinits(ngcx * ngcy, 10.0)
    a == "below18" && return [-8.0, 25.0, 30.0, density..., geo...]
    a == "18-25" && return [-7.0, 18.0, 20.0, density..., geo...]
    a == "25-30" && return [-7.0, 18.0, 20.0, density..., geo...]
    a == "30-50" && return [-7.5, 23.0, 30.0, density..., geo...]
    a == "50-65" && return [-7.5, 23.0, 30.0, density..., geo...]
    a == "above65" && return [-7.5, 23.0, 30.0, density..., geo...]
end

function bound(a, ndc, ngcx, ngcy)
    lbdensity = fill(-100.0, ndc)
    lbgeo = fill(-100.0, ngcx * ngcy)
    ubdensity = fill(100.0, ndc)
    ubgeo = fill(100.0, ngcx * ngcy)

    pastelb(c) = vcat(c, lbdensity..., lbgeo...)
    pasteub(c) = vcat(c, ubdensity..., ubgeo...)

    a == "below18" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 50.0])
    a == "18-25" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 30.0])
    a == "25-30" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 30.0])
    a == "30-50" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
    a == "50-65" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
    a == "above65" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
end
