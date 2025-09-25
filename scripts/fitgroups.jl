using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using Distributions, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using IterTools, Mooncake, GeoStats, GeoIO, CairoMakie, Suppressor

include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/analyze.jl")
include("../src/georbf.jl")
include("../src/densityrbf.jl")
include("../src/results.jl")
include("../src/diagplots.jl")
include("../src/model.jl")
include("../src/modelutils.jl")

function fit_years(a, ndc, ngcx, p = 1.0)
    ## 2003 has data issues
    allyears = vcat(2000:2002, 2004:2017)
##    allyears = [2010]
    if isfile("output/optim$a")
        results = deserialize("output/optim$a")
        fittedyears = [r.mdl.data.year for r in results]
        years = setdiff(allyears, fittedyears)
        println("output/optim$a already exists, fitting only $(years) years")
    else
        results = EstimationResult[]
        years = allyears
    end

    Threads.@threads for y in years
        println("Starting $a in $y")
        if a == "30-50"
            out = fit30to50(y, p, ndc, ngcx)
        else
            out = fitage(a, y, p, ndc, ngcx)
        end
        println("$a in $y done")
        push!(results, out)
        serialize("output/optim$(a)", results)
    end
    return results
end

function fitage(a, y, p, ndc, ngcx)
    mdl = makemodel(a, y, p, ndc, ngcx)
    inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    return @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
end

function fit30to50(y, p, ndc, ngcx)
    ## 30-50 is harder to fit so we use MLEs from fit of half the
    ## data, which works better. We maybe could use better inits and
    ## fit full data right away, but we likely would need very
    ## specific inits / weights for the RBFs
    mdl = makemodel("30-50", y, p / 2, ndc, ngcx)
    inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    out = @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));

    mdl = makemodel("30-50", y, p, ndc, ngcx)
    ## inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    inits = vec(out.chn[1, 1 : end - 4, 1].value.data)
    return  @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
end

function makemodel(a, y, p, ndc, ngcx)
    mdl = baseflow(
        load_data(a, y, p, "../data/"; only_positive = true, seed = 1234),
        normalize = false, ndc = ndc, ngcx = ngcx, kgeo = 2.0
    );
    return mdl
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
for a in ages
    fit_years(a, 16, 5, .1)
end
