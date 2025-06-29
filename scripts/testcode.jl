using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
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

out1 = testrun(norm, data, 1, 1, false);
## out2 = testrun(norm, data, 1, 15, false);
out3 = testrun(norm, data, 15, 1, false);
## out4 = testrun(norm, data, 28, 1, false);

## Let's see what normalization does for us
out4 = testrun(norm, data, 1, 1, true);
out5 = testrun(norm, data, 15, 1, true);
## out6 = testrun(norm, data, 28, 1, true, 20);


inits = Float64.(collect(out.out[1: (end-3)]))
out = @time estimate(norm, data, ndc = 28, ngc = 28, normalize = true,
                     inits = inits);

out.out
out.plt[5]
out.plt[6]

AD = ADTypes.AutoForwardDiff()
mdl = norm(data; ndc = 28, ngc = 28, normalize = true);
Random.seed!(123)
mles = Turing.maximum_likelihood(mdl.mdl; lb = mdl.lb, ub = mdl.ub, adtype = AD,
                                 reltol = 1e-1, maxiters = 5, show_trace = true,
                                 extended_trace = true)

chn = Turing.sample(mdl.mdl, NUTS(), 5, progress = true)
