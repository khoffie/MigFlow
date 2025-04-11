using Turing, Distributions, CSV, DataFrames, StatsBase, Random
using Plots, StatsPlots, LaTeXStrings, KernelDensity, CategoricalArrays
using NamedArrays

include("src/diagplots.jl")
include("src/othermodels.jl")
include("src/estimation.jl")
include("src/utils.jl")
include("src/diag.jl")
include("../src/loadgermdata.jl")
include("src/fullmodel.jl")

districts = CSV.read("../data/districts.csv", DataFrame)
addlrd(districts)
df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
joinlrd(df, districts)
df = year(age(df, "30-50"), 2016)
df = sample_flows(df, .1)

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
plts = [distplots(age(df, a)) for a in ages]
plot(plts..., size = (1000, 700))

function estiml(flows, districts, ds)
    mdat = gen_mdat(flows, districts; distscale = ds,
                    ndc = 1, ngc = 1)
    l = estimate(distonly, mdat).out["l"]
    return l
end

# dscales = Float64.([1 : 1 : 10; 11 : 10 : 201; 202 : 50 : 803])
# ls = [estiml(df, districts, ds) for ds in dscales]
df2 = CSV.read("insensitive.csv", DataFrame)
df2.l = df2.l ./ 100


norm(x) = (x .- minimum(x)) ./ (maximum(x) .- minimum(x))
df2.dnorm = norm(df2.d)
frac(x, b, c) = 1 ./ (1 .+ (b ./ x.^c))

@model function insens(l, d)
    b ~ Gamma(1, 1.0)
    c ~ Gamma(3, 1.0)
    s ~ Gamma(50, 10)
    ##p  = 1 ./ (1 .+ (b ./ (d)).^c)
    p = frac(d, b, c)
    c = 1e-5
    p = clamp.(p, c, 1-c)
    l ~ arraydist(Beta.(p .* s, (1 .- p) .* s))
    return p
end

chn = Turing.sample(insens(df2.l, df2.dnorm), NUTS(), 1000)
chn
preds = generated_quantities(insens(df2.l, df2.dnorm),
                             chn.value[end, [:b, :c, :s]],
                             ["b", "c", "s"])
b = chn.value[end, :b]
c = chn.value[end, :c]

plot(df2.dnorm, df2.l, label = "'Observed' fraction of 1/r desirability")
plot!(df2.dnorm, preds, label = "Predicted")
plot!(df2.dnorm, frac(df2.dnorm, b, c))
plot(Gamma(50, 10))
chn

dsample = StatsBase.wsample(df.dist, Weights(df.flows), 10^3)
kd = kde(dsample, Normal(0, 10))
m = maximum(kd.density)

plot(kd.x, kd.density ./ m, label = "flows")
plot!(kd.x, frac(norm(kd.x), b, c), label = "fraction 1/r")

density(frac(norm(dsample), b, c))
plot(df2.dnorm, frac(df2.dnorm, b, c))
density(frac(df2.dnorm, b, c))
