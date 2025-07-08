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

rbfinits(N, σ, t = 90.0) = clamp.(rand(MvNormal(zeros(N), σ^2 *I(N))), -t, t)

function fit_years(a, p)
    years = vcat(2000:2002, 2004:2017)
    results = NamedTuple[]
    inits = [-7.5, 18.0, 20.0, rbfinits(9, 40.0)..., rbfinits(8, 10.0)...]
    for y in years
        mdl = norm(load_data(a, y, p, "../data/";
                             only_positive = true,
                             seed = 1234, opf = false),
                   normalize = false, ndc = 9, ngcx = 2);
        out = @time estimate(mdl, optim_kwargs = (; show_trace = false,
                                                  inits = inits));
        println("$y done")
        push!(results, out)
    end
    serialize("output/optim$(a)", results)
    return results
end

ages = ["below18", "25-30", "30-50", "50-65", "above65"]
for a in ages
    fit_years(a, 1.0)
end

# post = postprocess(results, vcat(2000:2002, 2004:2017))
# post.df
# idx = [1, 2, 3, 15, 16, 17]
# plot(post.pls[1:3]...)
# plot(post.pls[4][idx]..., size = (1200, 900))
# plot(post.pls[5][idx]..., size = (1200, 900))
