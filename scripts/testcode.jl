using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, LaTeXStrings, Loess, ADTypes, KernelDensity,
    IterTools, GeoStats, GeoIO, CairoMakie, DynamicPPL, Serialization, StatProfilerHTML

include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/analyze.jl")
include("../src/georbf.jl")
include("../src/densityrbf.jl")
include("../src/results.jl")
include("../src/modelutils.jl")
include("../src/plotutils.jl")
include("../src/utils.jl")
include("../src/seplot.jl")

include("../models/baseflow.jl")
include("../models/baseflownormalized.jl")
include("../models/fundamental.jl")
include("../models/gravity.jl")

shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp");

results, ages = readresults(["fundamental", "norm", "gravity"]);
results, ages = readresults(["baseflow"]);
results.baseflow.age_18to25.y2010.mdl.meta
mdl = baseflownormalized(load_data("18-25", 2017, 1.0, "../data/"; only_positive = true);
                         ndc = 16, ngcx = 5)
# mdl = fundamental(load_data("18-25", 2010, 1.0, "../data/"; only_positive = true))

out = @time estimate(mdl; optim_kwargs = (reltol = 1e-2, ))
ana = analyze(out)
ana.fig
m, p = plotdtf(out)
p
m, p = plotgeo(out, shp, st)
p

out.mdl.meta
# ana = analyze(out)
# m, p = plotdtf(out)
# m, p = plotgeo(out, shp, st)
# p
# ana.fig
out2 = @time estimate(mdl; optim_kwargs = (; maxtime = 100, ))
out2.mdl.meta
# results, ages = readresults(["baseflownorm"]);
# r = results.baseflownormalized.agealized_below18.y2010
# m, p = plotdtf(r)
# m, p = plotgeo(r, shp, st)
# p


# struct TruncatedPoisson{T<:Real} <: DiscreteUnivariateDistribution λ::T end

# function Distributions.logpdf(d::TruncatedPoisson, k::Real)
#     if k < 1
#         return -Inf
#     else
#         λ = d.λ
#         ## return k * log(λ) - λ - logfactorial(k) - log1mexp(-λ)
#         return logpdf(Poisson(λ), k) - log(1 - exp(-λ))
#     end
# end

# function Distributions.rand(rng::AbstractRNG, d::TruncatedPoisson)
#     k = 0
#     while k < 1
#         k = rand(rng, Poisson(d.λ))
#     end
#     return k
# end


# lvals = 0.01:0.02:1.0

# fig = Figure(size = (400, 400));
# ax = Axis(fig[1, 1])

# plot!(ax, lvals, Distributions.logpdf.(TruncatedPoisson.(lvals), 1.0))
# plot!(ax, lvals, Distributions.logpdf.(Poisson.(lvals), 1.0))
# fig

# using DifferentiationInterface

# f(x) = Distributions.logpdf(TruncatedPoisson(x), 1)
# backend = AutoMooncake(; config=nothing)
# prep = prepare_gradient(f, backend, .1)
# gradient(f, prep, backend, 10.0)


# x = 0
# x / (1 - exp(-x))
