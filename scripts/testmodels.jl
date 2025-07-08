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
plotcoef(df, c) = (plot(df.year, df[!, c], title = c); scatter!(df.year, df[!, c]))

function postprocess(results, loopvec)
    params = results[1].name_map.parameters
    df = DataFrame([Symbol(p) => [results[i].chn[p].data[1]
                                  for i in eachindex(results)] for p in params])
    df.group = loopvec
    pgamma = plotcoef(df, :γ_raw)
    palpha = plotcoef(df, :α_raw)
    pphi = plotcoef(df, :ϕ_raw)
    pdelta = plotcoef(df, :δ_raw)
    pdens = [results[i].plts[5] for i in legitfits]
    pgeo = [results[i].plts[6] for i in legitfits]
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

chn1 = results[1].chn
data = load_data("18-25", years[1], p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);
mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);

vals = chn1.value.data[1 :  end - 1] .* 20
params = chn1.name_map.parameters[1 : end - 1]
inits = NamedTuple{Tuple(params)}(Tuple(vals))

sam = Turing.sample(mdl.mdl, Turing.Inference.HMC(.01, 200), 1; init_theta = [inits])
