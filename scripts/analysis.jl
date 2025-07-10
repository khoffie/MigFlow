using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using BenchmarkTools, IterTools, StatProfilerHTML, ReverseDiff
using Suppressor
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

function reorder(results)
    years = [Int(r.chn[:year].data[1]) for r in results]
    return results[sortperm(years)]
end

function coefdf(results)
    params = results[1].chn.name_map.parameters
    df = DataFrame([Symbol(p) => [r.chn[p].data[1] for r in results] for p in params])
    df.group = [Int(r.chn[:year].data[1]) for r in results]
    first = ["α_raw", "γ_raw", "ϕ_raw", "lp", "group"]
    last = setdiff(names(df), first)
    select!(df, vcat(first, last))
    return sort!(df, :group)
end

plotcoef(df, c) = (plot(df.group, df[!, c], title = c); scatter!(df.group, df[!, c]))

function resplots(results, coefdf)
    pgamma = plotcoef(coefdf, :γ_raw)
    palpha = plotcoef(coefdf, :α_raw)
    pphi = plotcoef(coefdf, :ϕ_raw)
    # pdelta = plotcoef(dfp, :δ_raw)
    pdens = [r.plts[5] for r in results]
    pgeo = [r.plts[6] for r in results]
    return [palpha, pgamma, pphi, pdens, pgeo]
end

results = reorder(deserialize("./output/optim50-65"));
df = coefdf(results)
plts = resplots(results, df)
idx = [1, 2, 3, 15, 16, 17]
plot(plts[1:3]...)
plot(plts[4][idx]..., size = (1200, 900))
plot(plts[5][idx]..., size = (1200, 900))
