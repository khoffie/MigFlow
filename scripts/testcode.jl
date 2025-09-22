using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots, Distributions
using CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using IterTools, Mooncake, Revise, GeoStats, GeoIO, CairoMakie
using StatsBase: coeftable
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
includet("../src/analyze.jl")
includet("../src/analyzegeo.jl")
includet("../src/analyzedensity.jl")
includet("../src/analyzeresults.jl")
includet("../src/diagplots.jl")
include("../src/model.jl")
include("../src/model_helpers.jl")
include("../src/plotutils.jl")

shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp")

mdl = baseflow(
    load_data(
        "30-50", # age group
        2014, # year
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

inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
@time out = estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
## diagnostic plots
post = analyze(out)
## heatmap of density transition function
m, pdtf = plotdtf(out)
## Map of Germany showing locational asymmetries
geo, pgeo = plotgeo(out, shp, st)

savefig(post.plts[end], "../docs/check.pdf")
save("../docs/pdtf.pdf", pdtf)
save("../docs/pgeo.pdf", pgeo)
