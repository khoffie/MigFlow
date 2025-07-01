using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using Enzyme

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
    scale = rbfscale(cx, cy, 2.0)
    mat = [interpolant(rbf, Rfrom, Rto, coefs, cx, cy, scale)
           for Rfrom in vals, Rto in vals]';
    mat = mat .- mean(mat)
    display(heatmap(mat))
    return mat
end

p = 1.0
p = 0.1
data = load_data("below18", 2017, p, "../data/"; only_positive = true,
                 seed = 1234, opf = false);

out = @time estimate(norm, data;
                     model_kwargs = (; ndc = 16, ngc = 1, normalize = false),
                     optim_kwargs = (; ad = AutoEnzyme()));
out.out
out.plt[6]


AD = ADTypes.AutoForwardDiff()
AD = AutoEnzyme()
mdl = norm(data; W = 16, densscale = 2.0, ndc = 1, ngc = 1, normalize = false);
Random.seed!(123)
chn = Turing.sample(mdl.mdl, NUTS(100, .5), 1, progress = true)
denscoefs = extract_coefs(chn[end, :, 1], "ζ")
geocoefs = extract_coefs(chn[end, :, 1], "η")

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

W = 16
xcoord     = scale_to_unit(districts.xcoord)
ycoord     = scale_to_unit(districts.ycoord)
xmin, xmax = extrema(xcoord)
ymin, ymax = extrema(ycoord)
cxgeo      = [range(xmin, xmax, Int(sqrt(W)));]
cygeo      = [range(ymin, ymax, Int(sqrt(W)));]
geoscale   = rbfscale(cxgeo, cygeo, 2.0)

val = [interpolant(rbf, xcoord[i], ycoord[i], geocoefs, cxgeo, cygeo, geoscale)
       for i in eachindex(xcoord)]
val = val .- mean(val)
dfgeo = DataFrame(; xcoord, ycoord, val)

ratio = (ymax - ymin) / (xmax - xmin)
width = 600
p = scatter(dfgeo.xcoord, dfgeo.ycoord,
            marker_z = dfgeo.val,
                markersize = 10, size = (600, width * ratio),
            label = "")
p
