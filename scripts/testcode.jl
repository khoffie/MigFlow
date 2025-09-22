using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots, Distributions
using CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using IterTools, Mooncake, Revise, GeoStats, GeoIO, CairoMakie

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
@time out2 = estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));
inits = vec(out2.chn[1, 1 : end - 4, 1].value.data)
@time out2 = estimate(mdl, optim_kwargs = (; show_trace = false, inits = inits));

df = CSV.read("../data/FlowDataGermans.csv", DataFrame)

setdiff(lc(year(df, 2000).todist), lc(year(df, 2000).fromdist))

lvls = unique(sort(df.todist))

levelcode.(categorical(df.fromdist, levels = lvls))
setdiff(lc(year(df, 2000).todist), levelcode.(categorical(year(df, 2000).fromdist, levels = lvls)))
