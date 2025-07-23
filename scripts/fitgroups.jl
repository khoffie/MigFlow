using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using BenchmarkTools, IterTools, StatProfilerHTML, ReverseDiff, Distributions
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
include("../src/model_helpers.jl")

function fit_years(a, ndc, ngcx)
    ## 2003 has data issues
    allyears = vcat(2000:2002, 2004:2017)
    if isfile("output/optim$a")
        results = deserialize("output/optim$a")
        fittedyears = [Int(r.chn[:year].data[1]) for r in results]
        years = setdiff(allyears, fittedyears)
        println("output/optim$a already exists, fitting only $(years) years")
    else
        results = NamedTuple[]
        years = allyears
    end

    Threads.@threads for y in years
        println("Starting $a in $y")
        if a == "30-50"
            out = fit30to50(y, ndc, ngxx)
        else
            out = fitage(a, y, ndc, ngcx)
        end
        println("$a in $y done")
        push!(results, out)
        serialize("output/optim$(a)", results)
    end
    return results
end

function fitage(a, y, ndc, ngcx)
    mdl = makemodel(a, y, ndc, ngcx)
    inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    return @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
end

function fit30to50(y, ndc, ngcx)
    ## 30-50 is harder to fit so we use MLEs from fit of half the
    ## data, which works better. We maybe could use better inits and
    ## fit full data right away, but we likely would need very
    ## specific inits / weights for the RBFs
    mdl = makemodel("30-50", y, .5, ndc, ngcx)
    inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    out = @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));

    mdl = makemodel("30-50", y, 1.0, ndc, ngcx)
    ## inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    inits = out.chn[1, 1 : end - 4, 1].value.data
    return  @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = reshape(inits, 20)));
end

function makemodel(a, y, p, ndc, ngcx)
    mdl = norm(load_data(a, y, p, "../data/";
                         only_positive = true,
                         seed = 1234, opf = false),
               normalize = false, ndc = ndc, ngcx = ngcx);
    return mdl
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
for a in ages
    fit_years(a, 9, 5)
end


function savepredictions(files)
    for f in files
        results = deserialize(f);
        res = NamedTuple[]
        for r in results
            p = invokelatest(Turing.returned, r.mdl.mdl, r.chn)[1]
            push!(res, (mdl = r.mdl, chn = r.chn, prd = p))
        end
        serialize(f, res)
    end
end

savepredictions(readdir("output"; join = true))

