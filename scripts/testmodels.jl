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

        data = load_data(a, y, p, "../data/"; only_positive = true,
                         seed = 1234, opf = false);
        mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);
        out = @time estimate(mdl, optim_kwargs = (; show_trace = false,
                                                  inits = inits));
        println("$y done")
        push!(results, out)
    end
    serialize("output/optim$(a)", results)
    return results
end

results = fit_years("18-25", 1.0)



out = postprocess(results, vcat(2000:2002, 2004:2017), -200e3);
out = postprocess(results, vcat(2000:2002, 2004:2017));
out.df
plot(results[3].plts[1:4]...)
out
results[2].plts[5]

plot(out.pls[1:3]...)
plot(out.pls[4][1:6]...)
plot(out.pls[4][7:13]...)

plot(out.pls[5][1:6]...)
plot(out.pls[5][1:6]...)
out.df[1, 6: 14]

ws = [coefmat(collect(out.df[i, 6:14])) for i in 1 : nrow(out.df)]
ws[1]
ws[2]
names(out.df)
out.df[12, 6:14]
out.df[12, 15:22]

ws = [coefmat(collect(out.df[i, 15 : 22])) for i in 1 : nrow(out.df)]
ws[8]

p = 1.0
p = 0.5
data = load_data("18-25", 2017, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);
mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);


out = estimate(mdl, optim_kwargs = (; show_trace = false,
                                    inits = inits));
plot(out.plts[1:4]...)
res1 = out
res1.chn



mles = try_optim_once(mdl.mdl, mdl.lb, mdl.ub)
mles.lp

runoptim(mdl.mdl, mdl.lb, mdl.ub; lp_thresh = -1)

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

out = postprocess(results, seeds, -200e3)
out
plot(results[5].plts[1:4]...)
plot(out.pls[1:3]...)

plot(out.pls[4][1:6]...)
plot(out.pls[5][1:3]..., layout = (3, 1))


p = .1
data = load_data("18-25", 2017, p, "../data/"; only_positive = true,
                 seed = seeds[1], opf = false);
mdl = norm(data, normalize = false, ndc = 9, ngcx = 2);
out = @time estimate(mdl, optim_kwargs = (; lp_thresh = -5000))


mdl
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
