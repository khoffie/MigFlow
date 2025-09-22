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

mdl = norm(load_data("30-50", 2014, 0.1, "../data/";
                     only_positive = true,
                     seed = 1234, opf = false),
           normalize = false, ndc = 16, ngcx = 5);
inits = initialize(mdl.data.age, mdl.mdl.args.ndc, mdl.mdl.args.ngcx, mdl.mdl.args.ngcy);
@time out = estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
