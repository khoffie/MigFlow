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

# out = deserialize("../output/anchored");
# df, net, figs = analyze(out)
# coefplot(out)
# m, pdtf = plotdtf(out)
# geo, pgeo = plotgeo(out, shp, st)
# extract_params(out)

# out.ses
# out2[17].ses
# out2.ses
# out2 = deserialize("./output/optim18-25");

# cdf, net, figs = analyze(out)
# coefplot(out)
# m, pdtf = plotdtf(out)
# geo, pgeo = plotgeo(out, shp, st)
# extract_params(out)
