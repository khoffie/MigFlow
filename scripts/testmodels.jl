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
include("../src/distonly.jl")

corround(x, y) = round(cor(x, y), digits = 2)
plotcoef(df, c) = (plot(df.idx, df[!, c], title = c); scatter!(df.idx, df[!, c]))

function postprocess(results, loopvec, lpthresh = -10e6)
    params = results[1].chn.name_map.parameters
    df = DataFrame([Symbol(p) => [results[i].chn[p].data[1]
                                  for i in eachindex(results)] for p in params])
    df.group = loopvec
    df.idx = 1:length(results)
    first = ["α_raw", "γ_raw", "ϕ_raw", "lp", "group"]
    last = setdiff(names(df), first)
    select!(df, vcat(first, last))

    idxsuccess = df.lp .> lpthresh
    dfp = df[idxsuccess, :]
    pgamma = plotcoef(dfp, :γ_raw)
    palpha = plotcoef(dfp, :α_raw)
    pphi = plotcoef(dfp, :ϕ_raw)
    # pdelta = plotcoef(dfp, :δ_raw)
    pdens = [r.plts[5] for r in results[idxsuccess]]
    pgeo = [r.plts[6] for r in results[idxsuccess]]
    pls = [palpha, pgamma, pphi, pdens, pgeo]
    return (; df, pls)
end



p = 1.0
p = 0.1
data = load_data("18-25", 2017, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);
mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);
out = estimate(mdl, optim_kwargs = (; show_trace = true));

years = vcat(2000:2002, 2004:2017)
results = NamedTuple[]
for y in years
    p = 1.0
    data = load_data("18-25", y, p, "../data/"; only_positive = true,
                     seed = 1234, opf = false);
    mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);
    out = estimate(mdl, optim_kwargs = (; show_trace = false));
    println("$y done")
    push!(results, out)
end

serialize("output/optim18-25", results)

seeds = rand(1:1000, 10)

results = NamedTuple[]
for s in seeds
    p = .5
    data = load_data("18-25", 2017, p, "../data/"; only_positive = true,
                     seed = s, opf = false);
    mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);
    out = @time estimate(mdl, optim_kwargs = (; show_trace = false));
    println("$s done")
    push!(results, out)
end

out = postprocess(results, seeds, -50e3)
plot(out.pls[1:3]...)
plot(out.pls[4][1:3]...)
plot(out.pls[5][1:3]...)


p = .1
data = load_data("18-25", 2017, p, "../data/"; only_positive = true,
                 seed = seeds[1], opf = false);
mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);
out = estimate(mdl)

out2 = estimate(mdl, optim_kwargs = (; reltol = 1e-10))
out2.chn
out.chn
plot(out.plts[1:4]...)
plot(out.plts[6])

vals = chn1.value.data[1 :  end - 1]
params = chn1.name_map.parameters[1 : end - 1]
inits = NamedTuple{Tuple(params)}(Tuple(vals))

sam = Turing.sample(mdl.mdl, Turing.Inference.HMC(.01, 200), 100; initial_params = vals)
sam[:lp]

chn1[:lp]
sam = Turing.sample(mdl.mdl, NUTS(), 1; init_params = [inits])
