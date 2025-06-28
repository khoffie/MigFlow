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

function cornet(out)
    p1 = scatter(out.net.nmr, out.net.nmrp)
    p2 = scatter(out.net.net, out.net.netp)
    display(plot(p1, p2))
    println("Cor nmr: $(cor(out.net.nmr, out.net.nmrp))")
    println("Cor net: $(cor(out.net.net, out.net.netp))")
end

p = 1.0
p = 0.1
data = load_data("18-25", 2016, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);

function testrun(norm, data, ndc, ngc, normalize, maxmin = 20)
    Random.seed!(123)
    out = @time estimate(norm, data;
                         model_kwargs = (; ndc = ndc, ngc = ngc,
                                         normalize = normalize),
                         optim_kwargs = (; show_trace = false,
                                     maxtime = maxmin * 60));
    f = "out_ndc$(ndc)_normalize$(normalize)"
    p = "/home/konstantin/code/scripts/output/"
    serialize(joinpath(p, f), out)
    println("$f saved")
    return out
end

out1 = testrun(norm, data, 1, 1, false);
out2 = testrun(norm, data, 15, 1, false);
out3 = testrun(norm, data, 28, 1, false);

out3.out
out3.plt[5]
sort(out3.net[:, [:nmr, :nmrp]], :nmr)

### 28 chebies is much better with both correlations, however,
### correlation between net and predicted net is still below zero!
## All with 10% sample of data
# julia> cornet(out)
# Cor nmr: 0.29709508420026887
# Cor net: -0.6645164962014762

# julia> cornet(out2)
# Cor nmr: 0.42324015023232686
# Cor net: -0.1310660226079404


## Let's see what normalization does for us
out4 = testrun(norm, data, 1, 1, true);
out5 = testrun(norm, data, 15, 1, true);
out6 = testrun(norm, data, 28, 1, true, 20);




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
