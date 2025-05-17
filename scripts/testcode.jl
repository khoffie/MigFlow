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

data = load_data("30-50", 2017, 1.0, "../data/"; positive =  false, full = true);

from = sort(unique(data.dffull.fromdist))
to = sort(unique(data.dffull.todist))
df = DataFrame(; data.dffull.fromdist, data.dffull.todist, data.dffull.dist)
distmat = reshape(df.dist, 400, 401)' ./ 100.0
P = log.(data.districts.pop ./ 153000)
to = levelcode.(categorical(data.dffull.todist))

desmat = [exp(desirability(P[i], distmat[i, j], 0.15, 0.14, 4)) for i in 1:401, j in 1:400]
denom = [sum(desmat[i,j] for j in 1:400) for i in 1:401]




[P[i] + log(distmat[i, j]^-4) for i in 1:401, j in 1:400]


distmat[1, :]
distmat[1, :] .== df[df.fromdist .== 1001, :].dist


desirability(P, D, ϕ, δ, γ) = P + log((ϕ + (1 - ϕ) / (D + δ)^γ))



out1 = @time estimate(distonly, data);
out1.out
out2 = @time estimate(distnorm, data);
out2.out

sum(logpdf.(Poisson.(out1.preds), data.df.flows))
sum(logpdf.(Poisson.(out2.preds), data.df.flows))


mdat = gen_mdat("30-50", 2017, 0.1, "../data"; type = "joint",
                distscale = 100.0, ndc = 28, ngc = 1, only_positive = true);
out1 = @time estimate(distonly, mdat);
out1.out



out1 = @time estimate(norm, mdat);
out1.out
out4 = @time estimate(choice, mdat; normalize = false);
out4.out


mdat = gen_mdat(data; type = "conditional", distscale = 100.0, ndc = 28, ngc = 1);
out1 = @time estimate(norm, mdat);
out2 = @time estimate(choice, mdat);
out3 = @time estimate(distonly, mdat);
out4 = @time estimate(choice, mdat; normalize = false);
out1.out
out2.out
out3.out
out4.out

plotfit(mdat.flows, out1.preds)
mdat.flows
smoother!(log.(out1.preds), log.(mdat.flows))


x = log.(out1.preds)
y = log.(mdat.flows)
mod = loess(x, y)
us = range(extrema(x)..., step = .01)
Loess.predict(mod, us)

plot(out1.plt[1], out3.plt[1])

sum(logpdf.(Poisson.(preds), flows))
out1.out
sum(logpdf.(Poisson.(out3.preds), mdat.flows))
out3.out
