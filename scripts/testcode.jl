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

di = year(CSV.read("../data/districts.csv", DataFrame), 2017)
di.lc = lc(di.distcode)
shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp");

results, ages = readresults(["baseflow"], "./output");
results, ages = readresults(["baseflow"], "./output_trunc_norm/");

r = results.baseflow.age_18to25.y2017;
r = results.baseflow_truncated_normalized.age_below18_true_true.y2017;
r.mdl.meta
rana = analyze(r);
rm, rp = plotdtf(r);
rm2, rp2 = plotgeo(r, shp, st);


mdl = baseflow(load_data("18-25", 2017, 0.1, "../data/"; only_positive = true);
                  ndc = 25, ngcx = 5, trunc = false, norm = false);
out = @time estimate(mdl);
inits = vcat(out.ses.coef[1], 1.0, out.ses.coef[2:end])
mdl = baseflow(load_data("18-25", 2017, 1.0, "../data/"; only_positive = true);
                  ndc = 25, ngcx = 5, trunc = true, norm = true);
out = @time estimate(mdl, optim_kwargs = (; inits = inits, maxiters = 1, maxtime = 600));

out.mdl.meta
ana = analyze(out);
m, p = plotdtf(out);
m2, p2 = plotgeo(out, shp, st);

p
rp
rp2
p2

fig = Figure((size = (800, 400)));
ax1 = Axis(fig[1, 1])
scatter!(ax1, df1.asymp, df1.asym)
smoother!(ax1, df1.asymp, df1.asym)

ax2 = Axis(fig[2, 1])
scatter!(ax2, df2.asymp, df2.asym)
smoother!(ax2, df2.asymp, df2.asym)

fig
sort(df1, :asym)
sort(df2, :asym)


ana.df
rana.df

scatter(df1.asymp, df2.asymp, color = df2.)

net = ana.net

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
