using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots, Distributions
using CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using IterTools, Mooncake, GeoStats, GeoIO, CairoMakie
using StatsBase: coeftable
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/analyze.jl")
include("../src/georbf.jl")
include("../src/densityrbf.jl")
include("../src/results.jl")
include("../src/model.jl")
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
        0.1, # Fraction of rows to use, e.g. 10%
        "../data/"; ## path where FlowDataGermans.csv and districts.csv
        ## are stored
        only_positive = true, # return only positive flows / drop zero
        # flows
        seed = 1234, # Random seed for reproducibility
        opf = false # depracated, ignore
    ),
    normalize = false, ## normalize desirabilities, ## currently only false supported
    ndc = 16, # number of radial basis centers for density transition function
    ngcx = 5 # number of radial basis centers for geographical
             # asymmetries in x direction. y direction is set
             # automatically
);

inits = initialize(mdl);
@time out = estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
df, net, pcheck = analyze(out); ## diagnostic plots
m, pdtf = plotdtf(out) ## heatmap of density transition function
geo, pgeo = plotgeo(out, shp, st) ## Map of locational asymmetries
pcoefs = coefplot(out)
save("../docs/pcheck.png", pcheck)
save("../docs/pdtf.png", pdtf)
save("../docs/pgeo.png", pgeo)
save("../docs/pcoefs.png", pcoefs)
