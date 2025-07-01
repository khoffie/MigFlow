using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
include("../src/rbf.jl")
## available models
include("../src/norm.jl")

function extract_coefs(chn, string)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end

corround(x, y) = round(cor(x, y), digits = 2)

function heat(coefs, districts)
    R = scale_to_unit(log.(districts.density ./ median(districts.density)))
    Rmin, Rmax = extrema(R)
    vals = range(Rmin, Rmax, 1000)

    cx = [range(Rmin, Rmax, 4);]
    cy = [range(Rmin, Rmax, 4);]

    mat = [interpolant(rbf, Rfrom, Rto, coefs, cx, cy)
           for Rfrom in vals, Rto in vals]';
    mat = mat .- mean(mat)
    display(heatmap(mat))
    return mat
end

p = 1.0
p = 0.1
data = load_data("below18", 2017, p, "../data/"; only_positive = true,
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
chn = Turing.sample(mdl.mdl, NUTS(100, .6), 10, progress = true)
coefs = extract_coefs(chn[end, :, 1], "ω")
plot(chn[:lp], label = "$(round(chn[:lp][end], digits = 2))")

preds = returned(mdl.mdl, chn[end])[1];
df = DataFrame(; data.df.fromdist, data.df.todist, data.df.flows, preds, data.df.dist);
plot(plotfit(df.flows, df.preds),
     title = "Cor: $(corround(log.(df.flows), log.(df.preds)))")

net = calc_net_df(df);
plotnet(net)

districts = year(CSV.read("../data/districts.csv", DataFrame), 2017)
mat = heat(coefs, districts)
heatmap(mat)
