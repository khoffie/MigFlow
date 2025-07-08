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
plotcoef(df, c) = (plot(df.group, df[!, c], title = c); scatter!(df.group, df[!, c]))

function postprocess(results, loopvec)
    params = results[1].chn.name_map.parameters
    df = DataFrame([Symbol(p) => [results[i].chn[p].data[1]
                                  for i in eachindex(results)] for p in params])
    df.group = loopvec
    pgamma = plotcoef(df, :γ_raw)
    palpha = plotcoef(df, :α_raw)
    pphi = plotcoef(df, :ϕ_raw)
    pdelta = plotcoef(df, :δ_raw)
    pdens = [r.plts[5] for r in results]
    pgeo = [r.plts[6] for r in results]
    pls = [palpha, pgamma, pdelta, pphi, pdens, pgeo]
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
    p = 0.1
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
    p = 0.2
    data = load_data("18-25", 2017, p, "../data/"; only_positive = true,
                     seed = s, opf = false);
    mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);
    out = @time estimate(mdl, optim_kwargs = (; show_trace = false));
    println("$s done")
    push!(results, out)
end

out = postprocess(results, seeds)

sort(out.df, :lp)[!, vcat(1:5, end - 1, end)]



chn1 = results[1].chn
data = load_data("18-25", 2017, p, "../data/"; only_positive = true,
                 seed = seeds[1], opf = false);
mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);

vals = chn1.value.data[1 :  end - 1]
params = chn1.name_map.parameters[1 : end - 1]
inits = NamedTuple{Tuple(params)}(Tuple(vals))

sam = Turing.sample(mdl.mdl, Turing.Inference.HMC(.01, 200), 100; initial_params = vals)
sam[:lp]
chn1[:lp]
sam = Turing.sample(mdl.mdl, NUTS(), 1; init_params = [inits])
