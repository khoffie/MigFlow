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
include("../src/diagchain.jl")
include("../src/diagrbf.jl")
include("../src/norm.jl")
include("../src/rbf.jl")

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
    ## ub alpha only makes sense for distscale = 100 and pop /
    ## median(pop). Otherwise base prob to migrate might be very different
    a == "below18" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 50.0])
    a == "18-25" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 30.0])
    a == "25-30" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 30.0])
    a == "30-50" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
    a == "50-65" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
    a == "above65" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
end

function fit_years(a, p)
    years = vcat(2000:2002, 2004:2017)
    results = NamedTuple[]

    Threads.@threads for y in years
        mdl = norm(load_data(a, 2014, p, "../data/";
                             only_positive = true,
                             seed = 1234, opf = false),
                   normalize = false, ndc = 9, ngcx = 2);
        inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
        out = @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));

        println("$y done")
        push!(results, out)
        serialize("output/optim$(a)", results)
    end
    return results
end

ages = ["below18", "25-30", "30-50", "50-65", "above65"]
for a in ages
    fit_years(a, 1.0)
end
