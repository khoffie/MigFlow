using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, LaTeXStrings, Loess, ADTypes, KernelDensity,
    IterTools, GeoStats, GeoIO, CairoMakie, DynamicPPL, Serialization, SpecialFunctions,
    LogExpFunctions, LinearAlgebra

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

include("./maplots.jl")

results, ages = readresults(["fundamental", "gravity", "fundamental_normalized"], "./output");

mdl = fundamental(load_data("25-30", 2017, 1.0, "../data/"; only_positive = true);
                  trunc = true, norm = false);
mdl = fundamental(load_data("25-30", 2017, 1.0, "../data/"; only_positive = true);
                  trunc = false, norm = false);
out = estimate(mdl)
ana = analyze(out)
ana.fig
df = ana.df


res = results.gravity.age_18to25.y2000;
res = results.gravity.age_30to50.y2017;
res = results.fundamental.age_18to25.y2000;
res = results.fundamental.age_30to50.y2017;

analyze(res).fig
plotmflows(res)
predhist(res, 2)
plotmdist(res)
