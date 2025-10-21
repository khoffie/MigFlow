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

mdl = normalized(load_data("50-65", 2017, 1.0, "../data/"; only_positive = true); trunc = false);
out = @time estimate(mdl)
ana = analyze(out)

mdl2 = normalized(load_data("50-65", 2017, 1.0, "../data/"; only_positive = true); trunc = true);
out2 = @time estimate(mdl2)
ana2 = analyze(out2)

mdl3 = fundamental(load_data("50-65", 2017, 1.0, "../data/"; only_positive = true));
out3 = @time estimate(mdl3)
ana3 = analyze(out3)

mdl4 = baseflownormalized(load_data("50-65", 2017, 1.0, "../data/"; only_positive = true);
                          ndc = 16, ngcx = 5, trunc = true);
coefs = r.ses.coef
inits = vcat(coefs[1], 1.0, coefs[2:end]) # need init for beta
out4 = @time estimate(mdl4; optim_kwargs = (; inits = inits))
ana4 = analyze(out4)

out.mdl.meta
out2.mdl.meta
out3.mdl.meta
ana.fig
ana2.fig
ana3.fig
