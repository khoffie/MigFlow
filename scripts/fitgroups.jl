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

rbfinits(N, σ, t = 90.0) = clamp.(rand(MvNormal(zeros(N), σ^2 *I(N))), -t, t)

function initialize(a, ndc, ngcx, ngcy)
    density = rbfinits(ndc, 40.0)
    geo = rbfinits(ngcx * ngcy, 10.0)
    a == "below18" && return [-8.0, 25.0, 30.0, density..., geo...]
    a == "18-25" && return [-7.0, 18.0, 20.0, density..., geo...]
    a == "25-30" && return [-7.0, 18.0, 20.0, density..., geo...]
    a == "30-50" && return [-6.5, 20.0, 18.0, density..., geo...]
    a == "50-65" && return [-7.5, 23.0, 30.0, density..., geo...]
    a == "above65" && return [-7.5, 23.0, 30.0, density..., geo...]
end

function fit_years(a)
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
            out = fit30to50(y)
        else
            out = fitage(a, y)
        end
        println("$a in $y done")
        push!(results, out)
        serialize("output/optim$(a)", results)
    end
    return results
end

function fitage(a, y)
    mdl = norm(load_data(a, y, 1.0, "../data/";
                         only_positive = true,
                         seed = 1234, opf = false),
               normalize = false, ndc = 9, ngcx = 2);
    inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    return @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
end

function fit30to50(y)
    ## 30-50 is harder to fit so we use MLEs from fit of half the
    ## data, which works better. We maybe could use better inits and
    ## fit full data right away, but we likely would need very
    ## specific inits / weights for the RBFs
    mdl = norm(load_data("30-50", y, 0.5, "../data/";
                         only_positive = true,
                         seed = 1234, opf = false),
               normalize = false, ndc = 9, ngcx = 2);
    inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    out = @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));

    mdl = norm(load_data("30-50", y, 1.0, "../data/";
                         only_positive = true,
                         seed = 1234, opf = false),
               normalize = false, ndc = 9, ngcx = 2);
    ## inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
    inits = out.chn[1, 1 : end - 4, 1].value.data
    return  @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = reshape(inits, 20)));
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
for a in ages
    fit_years(a)
end
