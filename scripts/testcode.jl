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
include("../models/TruncatedPoisson.jl")


mdl = baseflow(load_data("25-30", 2017, .1, "../data/"; only_positive = true);
                  ndc = 25, ngcx = 5, trunc = false, norm = false);
out = @time estimate(mdl; optim_kwargs = (; maxtime = 100))
## inits = vcat(out.ses.coef[1], 1.0, out.ses.coef[2:end])
inits = out.ses.coef
mdl = baseflow(load_data("25-30", 2017, 1.0, "../data/"; only_positive = true);
                  ndc = 25, ngcx = 5, trunc = true, norm = false);

out = @time estimate(mdl; optim_kwargs = (; maxtime = 200, initial_params = inits))
out = @time estimate(mdl; optim_kwargs = (; maxtime = 200))

out = @time estimate(mdl; optim_kwargs = (; maxtime = 600, initial_params = out.ses.coef))
out.mdl.meta
out.ses
mdl.ub
a = analyze(out)
a.fig
