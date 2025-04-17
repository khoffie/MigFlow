using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings
using BenchmarkTools, ADTypes, ReverseDiff

include("src/utils.jl")
include("src/othermodels.jl")
include("src/estimation.jl")
include("src/loadgermdata.jl")
include("src/diag.jl")
include("src/diagplots.jl")
include("src/fullmodel.jl")

sizemb(x) =  Base.summarysize(x) / 1024^2

function load_data(a, y, p)
    di = CSV.read("data/districts.csv", DataFrame)
    di = add_lrd(di)
    df = CSV.read("data/FlowDataGermans.csv", DataFrame)
    df = year(age(pos(df), a), y)
    df = sample_flows(df, p)
    df = joinlrd(df, di)
    return (df = df, districts = di[di.year .== y, :])
end

function benchmark_model(data, ndc, ngc)
    b = Dict{Tuple{Int, Int}, BenchmarkTools.Trial}()
    for (i, j) in zip(ndc, ngc)
        mdat = gen_mdat(data; distscale = 100.0, ndc = i, ngc = j)
        println("Starting benchmark for $((i, j)) coefs")
        b[i, j] = @benchmark(estimate(full, $mdat), samples = 1)
        println("Time elabsed $(mean(b[i, j]))")
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

# ndc = [1, 4, 8, 16]
# ngc = [1, 1, 1, 1]
# b01 = benchmark_model(load_data("30-50", 2017, 0.1), ndc, ngc)
# plot_time(ncoefs, b01)

mdat = gen_mdat(load_data("30-50", 2017, 0.1);
                distscale = 100.0, ndc = 28, ngc = 28)
out = @time estimate(full, mdat);

out2 = @time estimate(distonly, mdat);
out2.out
out.out[1:5]
out.plt[3]
out.plt[4]
out.plt[5]
out.plt[6]


plot(1:36, out.out[6 : end-1], title = "Geocheby coefs, scaled as coef / 100")
scatter!(1:36, out.out[6 : end-1])


# out.plt[5]
# out.plt[6]

# data = load_data("30-50", 2017, 0.1)
# data
# Base.summarysize(data)
# df = CSV.read("data/FlowDataGermans.csv", DataFrame)
# Base.summarysize(data) / Base.summarysize(df)
# Base.summarysize(data)
