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
    a == "30-50" && return [-7.5, 23.0, 30.0, density..., geo...]
    a == "50-65" && return [-7.5, 23.0, 30.0, density..., geo...]
    a == "above65" && return [-7.5, 23.0, 30.0, density..., geo...]
end

function fit_years(a, p)
    allyears = vcat(2000:2002, 2004:2017)
    if isfile("output/optim$a")
        results = deserialize("output/optim$a")
        fittedyears = [Int(r.chn[:year].data[1]) for r in results]
        years = setdiff(allyears, fittedyears)
    else
        results = NamedTuple[]
        years = allyears
    end

    Threads.@threads for y in years
        println("Starting $a in $y")
        mdl = norm(load_data(a, y, p, "../data/";
                             only_positive = true,
                             seed = 1234, opf = false),
                   normalize = false, ndc = 9, ngcx = 2);
        inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
        out = @time estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));

        println("$a in $y done")
        push!(results, out)
        serialize("output/optim$(a)", results)
    end
    return results
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
for a in ages
    fit_years(a, 1.0)
end
