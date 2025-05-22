using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings
using ADTypes, ReverseDiff, Loess

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
include("../src/gen_mdat.jl")
## available models
include("../src/choiceset.jl")
include("../src/fullmodel.jl")
include("../src/fullmodel2.jl")
include("../src/othermodels.jl")
include("../src/norm.jl")

ratio(out) = exp(out.out["β_raw"] / exp(out.out["α_raw"]))
data = load_data("30-50", 2017, 1.0, "../data/"; only_positive = true,
                 seed = 1234);

out1 = @time estimate(norm, data; norm = "none", type = "conditional",
                      norm_type = "none");
n = "both"
out2 = @time estimate(norm, data; norm = n, type = "conditional",
                      norm_type = "sample");
out3 = @time estimate(norm, data; norm = n, type = "conditional",
                      norm_type = "full");

data = load_data("30-50", 2017, 0.1, "../data/"; only_positive = true);

ty = "joint"
out1 = @time estimate(norm, data; norm = "none", type = ty, norm_type = "none");
n = "both"
out2 = @time estimate(norm, data; norm = n, type = ty, norm_type = "sample");
out3 = @time estimate(norm, data; norm = n, type = ty, norm_type = "full");

lc(data.dffull.todist)


out1.out # conditional, no normalization
out2.out # conditional, normalization with subset
out3.out # conditional, normalization with full data

out1 = @time estimate(norm, data; norm = "origin", type = "conditional");
out1 = @time estimate(norm, data; norm = "destination", type = "conditional");
out1 = @time estimate(norm, data; norm = "both", type = "joint",
                      norm_type = "sample");



out1 = @time estimate(norm, data; norm = "none", type = "joint", norm_type = "none");
out2 = @time estimate(norm, data; norm = "both", type = "joint",
                      norm_type = "full");

out1.out # no normalization
out2.out # normalization


ratio(out1)
out1 = @time estimate(norm, data; norm = "origin", type = "joint");
out1 = @time estimate(norm, data; norm = "destination", type = "joint");

out1 = @time estimate(norm, data; norm = false, type = "joint");

D = data.dffull.dist
R = fradius.(data.districts.pop, data.districts.density)

D[D .< 30]

density(D)
density!(R)
plot(histogram(D), histogram(R))
exp(-1.3)
