using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, LaTeXStrings, Loess, ADTypes, KernelDensity,
    IterTools, GeoStats, GeoIO, CairoMakie, DynamicPPL, Serialization

include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/analyze.jl")
include("../src/georbf.jl")
include("../src/densityrbf.jl")
include("../src/results.jl")
include("../src/modelutils.jl")
include("../src/plotutils.jl")
include("../src/utils.jl")
include("../src/coefplot.jl")

include("../models/baseflow.jl")

shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp");

mdl = baseflow(
    load_data(
        "18-25", # age group
        2017, # year
        1.0, # Fraction of rows to use, e.g. 10%
        "../data/"; ## path of FlowDataGermans.csv and districts.csv
        only_positive = true, # use only positive flows / drop zero flows
        seed = 1234, # for reproducibility when sampling rows
    ),
    ndc = 16, # number of radial basis centers for density transition function
    ngcx = 5 # number of radial basis centers for geographical
             # asymmetries in x direction. y direction is set
             # automatically
);

serialize("./output/anchored", @time estimate(mdl))

out = deserialize("./output/anchored");
rename!(out.ses, [:coef => :coef_a, :se => :se_a])
out2 = deserialize("./output/optim18-25")[17];
df = leftjoin(out.ses, out2.ses, on = :name, makeunique = true)
df.diff = df.coef_a .- df.coef
fig = Figure(size = (1000, 400));
ax = Axis(fig[1, 1],
          xticks = (1:54, df.name),
          xgridvisible = false, ygridvisible = false,
          title = "Differences in estimates: anchored - non-anchored",
          ylabel = "Difference")
lines!(ax, 1:54, df.diff)
hlines!(ax, 0, color = :darkred)
fig
df
fig = Figure(size = (800, 400));
plotdtf(out, (-1, 1), fig, 1, 1, false)
plotdtf(out2, (-1, 1), fig, 1, 2)
prettytitle!(fig, "Left: Anchord, Right: Not anchored", -1)
fig

fig = Figure(size = (600, 400));
plotgeo(out, shp, st, (-.2, .6), fig, 1, 1, false)
plotgeo(out2, shp, st, (-.2, .6), fig, 1, 2)
prettytitle!(fig, "Left: Anchord, Right: Not anchored", -1)
fig
