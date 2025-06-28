using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
## available models
include("../src/norm.jl")

p = .1
p = 1.0
data = load_data("30-50", 2016, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);
Random.seed!(123)
out = @time estimate(norm, data, ndc = 28, ngc = 28, normalize = true);
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


using Dates


result = nothing
timelimit = 10.0  # seconds

task = @async begin
    result = Turing.maximum_likelihood(mdl.mdl;
        lb = mdl.lb, ub = mdl.ub,
        adtype = AD,
        reltol = 1e-6,
        maxiters = 10_000,
        callback = cb,
    )
end

timer = Timer(timelimit) do _
    Base.throwto(task, InterruptException())
end

try
    wait(task)
catch e
    @warn "Optimization interrupted after $timelimit seconds."
end
