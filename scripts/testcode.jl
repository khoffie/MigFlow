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
include("../models/baseflownormalized.jl")
include("../models/baseflowtruncated.jl")
include("../models/TruncatedPoisson.jl")

shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp");

results, ages = readresults(["baseflow"], "./output");
r = results.baseflow.age_50to65.y2017;

mdl = fundamental(load_data("50-65", 2017, 1.0, "../data/"; only_positive = true);
                  trunc = true, norm = true);
out = @time estimate(mdl);


mdl = baseflow(load_data("50-65", 2017, .1, "../data/"; only_positive = true);
                  ndc = 16, ngcx = 5, trunc = false, norm = false);
out = @time estimate(mdl);
