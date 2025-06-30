using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
include("../src/rbf.jl")
## available models
include("../src/norm.jl")


p = 1.0
p = 0.1
data = load_data("18-25", 2016, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);

function testrun(norm, data, ndc, ngc, normalize, maxmin = 20, save = false)
    Random.seed!(123)
    out = @time estimate(norm, data;
                         model_kwargs = (; ndc = ndc, ngc = ngc,
                                         normalize = normalize),
                         optim_kwargs = (; show_trace = false,
                                     maxtime = maxmin * 60));
    if save
        f = "out_ndc$(ndc)_normalize$(normalize)"
        p = "/home/konstantin/code/scripts/output/"
        serialize(joinpath(p, f), out)
        println("$f saved")
    end
    return out
end


## out2 = testrun(norm, data, 1, 15, false);
out3 = testrun(norm, data, 15, 1, false);
## out4 = testrun(norm, data, 28, 1, false);

## Let's see what normalization does for us
out4 = testrun(norm, data, 1, 1, true);
out5 = testrun(norm, data, 15, 1, true);
## out6 = testrun(norm, data, 28, 1, true, 20);


AD = ADTypes.AutoForwardDiff()
mdl = norm(data; W = 16, ndc = 1, ngc = 1, normalize = false);
Random.seed!(123)
chn = Turing.sample(mdl.mdl, NUTS(100, .6), 5, progress = true)
chn[:lp][end]
plot(chn[:lp])

DynamicPPL.DebugUtils.model_warntype(mdl.mdl)

out = @time estimate(norm, data; model_kwargs = (; W = 16, ndc = 1, ngc = 1,
                                                 normalize = false));
