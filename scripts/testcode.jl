using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, LaTeXStrings, Loess, ADTypes, KernelDensity,
    IterTools, GeoStats, GeoIO, CairoMakie, DynamicPPL, Serialization, SpecialFunctions,
    LogExpFunctions, StatProfilerHTML


include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/analyze.jl")
include("../src/georbf.jl")
include("../src/densityrbf.jl")
include("../src/results.jl")
include("../src/modelutils.jl")
include("../src/plotutils.jl")
include("../src/utils.jl")
include("../src/seplot.jl")
include("../src/plotcovcor.jl")

include("../models/fundamental.jl")
include("../models/baseflow.jl")
include("../models/gravity.jl")
include("../models/TruncatedPoisson.jl")


mdl = baseflow(load_data("18-25", 2017, 0.1, "../data/"; only_positive = true);
                  ndc = 25, ngcx = 5, trunc = false, norm = false);
@time maximum_a_posteriori(mdl.mdl; maxtime = 10)

out = @time estimate(mdl; optim_kwargs = (; maxtime = 100));
out.mdl.meta

net = vcat(ana.net, rana.net)
crange = extrema(vcat(net.asym, net.asymp))

function plotmap(ax, df, shp, st, col, crange = nothing)
    if isnothing(crange); crange = extrema(df[!, col]); end
    viz!(ax, shp.geometry, color = df[!, col], colorrange = crange);
    hideall!(ax)
    overlay_states(ax, st)
end

btn = "baseflow_truncated_normalized"
fig = Figure((size = (1200, 400)));
ax1 = Axis(fig[1, 1], aspect=DataAspect(), title = "normalized")
plotmap(ax1, selmodel(net, btn), shp, st, :asymp, crange);
ax2 = Axis(fig[1, 2], aspect=DataAspect(), title = "actual")
plotmap(ax2, selmodel(net, btn), shp, st, :asym, crange);
ax3 = Axis(fig[1, 3], aspect=DataAspect(), title = "non-normalized")
plotmap(ax3, selmodel(net, "baseflow"), shp, st, :asymp, crange);
Colorbar(fig[1, 4], limits = crange)
rowgap!(fig.layout, -50)
colgap!(fig.layout, -100)
fig

r.mdl.meta
out.mdl.meta
