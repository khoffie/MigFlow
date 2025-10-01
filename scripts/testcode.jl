using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, LaTeXStrings, Loess, ADTypes, KernelDensity,
    IterTools, GeoStats, GeoIO, CairoMakie, DynamicPPL

include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/analyze.jl")
include("../src/georbf.jl")
include("../src/densityrbf.jl")
include("../src/results.jl")
include("../src/baseflow.jl")
include("../src/fundamental.jl")
include("../src/gravity.jl")
include("../src/modelutils.jl")
include("../src/plotutils.jl")
include("../src/utils.jl")
include("../src/coefplot.jl")

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

mdl = fundamental(
    load_data(
        "18-25", # age group
        2017, # year
        1.0, # Fraction of rows to use, e.g. 10%
        "../data/"; ## path of FlowDataGermans.csv and districts.csv
        only_positive = true, # use only positive flows / drop zero flows
        seed = 1234, # for reproducibility when sampling rows
    )
);

@time out = estimate(mdl);
df, net, figs = analyze(out); ## diagnostic plots
figs
m, pdtf = plotdtf(out) ## heatmap of density transition function
geo, pgeo = plotgeo(out, shp, st) ## Map of locational asymmetries
pcoefs = coefplot(out)

mdl2 = gravity(
    load_data(
        "18-25", # age group
        2017, # year
        1.0, # Fraction of rows to use, e.g. 10%
        "../data/"; ## path of FlowDataGermans.csv and districts.csv
        only_positive = true, # use only positive flows / drop zero flows
        seed = 1234, # for reproducibility when sampling rows
    )
);

@time out2 = estimate(mdl2);
df2, net2, figs2 = analyze(out2); ## diagnostic plots
figs2
pcoefs2 = coefplot(out2)
pcoefs

mdl3 = norm(
    load_data(
        "18-25", # age group
        2017, # year
        1.0, # Fraction of rows to use, e.g. 10%
        "../data/"; ## path of FlowDataGermans.csv and districts.csv
        only_positive = true, # use only positive flows / drop zero flows
        seed = 1234, # for reproducibility when sampling rows
    )
);

@time out3 = estimate(mdl3);
df2, net2, figs2 = analyze(out2); ## diagnostic plots
figs2
pcoefs2 = coefplot(out2)
pcoefs
