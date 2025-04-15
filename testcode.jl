using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings
using BenchmarkTools

include("src/utils.jl")
include("src/othermodels.jl")
include("src/estimation.jl")
include("src/loadgermdata.jl")
include("src/diag.jl")
include("src/diagplots.jl")
include("src/fullmodel.jl")


function load_data(a, y, p)
    di = CSV.read("data/districts.csv", DataFrame)
    di = add_lrd(di)
    df = CSV.read("data/FlowDataGermans.csv", DataFrame)
    df = year(age(pos(df), a), y)
    df = sample_flows(df, p)
    df = joinlrd(df, di)
    return (df = df, districts = di[di.year .== y, :])
end

function benchmark_model(data, ncoefs)
    b = Dict{Int, BenchmarkTools.Trial}()
    for n in ncoefs
        mdat = gen_mdat(data; distscale = 100.0, ndc = n, ngc = 1)
        println("Starting benchmark for $n coefs")
        b[n] = @benchmark(estimate(full, $mdat), samples = 1)
        println("Time elabsed $(mean(b[n]))")
    end
    return b
end

function _plot_time(f, ncoefs, b)
    p = f(ncoefs, [mean(b[n]).time for n in ncoefs])
    scatter!(ncoefs, [mean(b[n]).time for n in ncoefs])
    return p
end

plot_time(ncoefs, b) = _plot_time(plot, ncoefs, b)
plot_time!(ncoefs, b) = _plot_time(plot!, ncoefs, b)

ncoefs = [1, 4, 8, 16, 25]
b01 = benchmark_model(load_data("30-50", 2017, 0.1), ncoefs)
plot_time(ncoefs, b01)
