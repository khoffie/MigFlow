using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings
using ADTypes, ReverseDiff, Loess, StatProfilerHTML

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
p = 1.0
n = "both"
ty = "joint"

out3f = @profilehtml estimate(norm,
                              load_data("30-50", 2017, p,
                                        "data/";
                                        only_positive = true,
                                        seed = 1234,
                                        opf = false);
                              norm = n, type = ty,
                              norm_type = "full");

# out1f = @time estimate(norm, data; norm = "none", type = ty, norm_type = "none");
# out2f = @time estimate(norm, data; norm = n, type = ty, norm_type = "sample");
