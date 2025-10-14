using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, LaTeXStrings, Loess, ADTypes, KernelDensity,
    IterTools, GeoStats, GeoIO, CairoMakie, DynamicPPL, Serialization, StatProfilerHTML

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

include("../models/baseflow.jl")
include("../models/baseflownormalized.jl")
include("../models/fundamental.jl")
include("../models/gravity.jl")

shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp");

results, ages = readresults(["fundamental", "norm", "gravity"]);

mdl = baseflownormalized(load_data("18-25", 2010, 1.0, "../data/"; only_positive = true);
                         ndc = 16, ngcx = 5)
## mdl = fundamental(load_data("18-25", 2010, 1.0, "../data/"; only_positive = true))

@profilehtml estimate(mdl)

# ana = analyze(out)
# m, p = plotdtf(out)
# m, p = plotgeo(out, shp, st)
# p
# ana.fig


# results, ages = readresults(["baseflownorm"]);
# r = results.baseflownormalized.agealized_below18.y2010
# m, p = plotdtf(r)
# m, p = plotgeo(r, shp, st)
# p
