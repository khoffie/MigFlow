using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using BenchmarkTools, IterTools

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
include("../src/norm.jl")
include("../src/rbf.jl")

function chaindiag(chain, mdl, data)
    chn[:lp][end]
    p1 = plot(chn[:lp], label = "$(round(chn[:lp][end], digits = 2))")

    denscoefs = extract_coefs(chn[end, :, 1], "ζ")
    mat, p2 = plotdensrbf(denscoefs, mdl.data)

    geocoefs = extract_coefs(chn[end, :, 1], "η")
    dfgeo, p3 = plotgeorbf(geocoefs, mdl.data)

    df = DataFrame(;
                   data.df.fromdist,
                   data.df.todist,
                   data.df.flows,
                   preds = returned(mdl.mdl, chn[end])[1],
                   data.df.dist);

    p4 = plot(plotfit(df.flows, df.preds),
              title = "Cor: $(corround(log.(df.flows), log.(df.preds)))")

    p5 = plotdist(df, :preds)

    net = calc_net_df(df);
    p6 = plot(plotnet(net), title = "cor: $(corround(net.nmr, net.nmrp))")
    plts = [p1, p2, p3, p4, p5, p6]
    return (; df, net, plts)
end

function extract_coefs(chn, string)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end

corround(x, y) = round(cor(x, y), digits = 2)

p = 1.0
p = 0.1
data = load_data("below18", 2017, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);

out = @time estimate(norm, data;
                     model_kwargs = (; ndc = 4, ngc = 4, normalize = false),
                     optim_kwargs = (; ad = AD));


AD = AutoEnzyme()
AD = ADTypes.AutoForwardDiff()
mdl = norm(data; ndc = 4, ngc = 4, normalize = false);
chn = Turing.sample(mdl.mdl, NUTS(10, .7), 1, progress = true)
out = chaindiag(chn, mdl, data)
## tree depth 4 or 5
## check how NUTS works
## rescale y axis, such that 1m = 1m
