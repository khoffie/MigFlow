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
include("../models/fundamental.jl")
include("../models/gravity.jl")

shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp");

mdl = baseflow(
    load_data(
        "25-30", # age group
        2002, # year
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

out = @time estimate(mdl);
df, net, quick, figs = analyze(out)
m, p1 = plotdtf(out)
g, p2 = plotgeo(out, shp, st)
figs
p1
p2

p3  = plotcovcor(out)

fig2 = Figure(size = (1000, 400));
ax = Axis(fig2[1, 1],
          xticks = (1:53, nms[2:end]),
          xgridvisible = false, ygridvisible = false,
          title = "correlations with alpha",
          xticklabelrotation = -π/4)
barplot!(ax, 1:53, cormat[1, 2 : end].array)
fig2
cormat
