using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, LaTeXStrings, Loess, ADTypes, KernelDensity,
    IterTools, GeoStats, GeoIO, CairoMakie, DynamicPPL, Serialization, SpecialFunctions,
    LogExpFunctions

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
include("../src/plotcovcor.jl")

include("../models/baseflow.jl")
include("../models/baseflowtruncated.jl")
include("../models/fundamental.jl")
include("../models/gravity.jl")
include("../models/TruncatedPoisson.jl")

shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp");

results, ages = readresults(["fundamental", "norm", "gravity", "baseflow"]);

mdl = baseflowtruncated(load_data("18-25", 2010, 1.0,
                                  "../data"; only_positive = true);
                        ndc = 16, ngcx = 5);


mdl1 = truncated(load_data("18-25", 2010, 1.0, "../data"; only_positive = true); type = 1);
mdl2 = truncated(load_data("18-25", 2010, 1.0, "../data"; only_positive = true); type = 2);
mdl3 = truncated(load_data("18-25", 2010, 1.0, "../data"; only_positive = true); type = 3);

out = @time estimate(mdl; optim_kwargs = (inits = r.ses.coef,));
