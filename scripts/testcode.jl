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

function plotgeo(coefs, districts)
    xcoord     = scale_to_unit(districts.xcoord)
    ycoord     = scale_to_unit(districts.ycoord)
    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    cxgeo      = [range(-0.6, 0.4, Int(sqrt(length(coefs))));]
    cygeo      = [range(-0.9, 0.6, Int(sqrt(length(coefs))));]
    geoscale   = rbfscale(cxgeo, cygeo, 1.0)
    geo = [interpolant(rbf, xcoord[i], ycoord[i], coefmat(coefs),
                       cxgeo, cygeo, geoscale) for i in eachindex(xcoord)];
    geo = geo .- mean(geo)
    dfgeo = DataFrame(; xcoord, ycoord, geo)
    ratio = (ymax - ymin) / (xmax - xmin)
    width = 600
    p = scatter(xcoord, ycoord,
                marker_z = dfgeo.geo,
                markersize = 10, size = (600, width * ratio),
                label = "")
    res = collect(product(cxgeo, cygeo))
    xvals = getindex.(res, 1)
    yvals = getindex.(res, 2)
    scatter!(vec(xvals), vec(yvals), color = :red)
    display(p)
    return df
end

function extract_coefs(chn, string)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end

corround(x, y) = round(cor(x, y), digits = 2)

function heat(coefs, k, districts)
    R = scale_to_unit(log.(districts.density ./ median(districts.density)))
    Rmin, Rmax = extrema(R)
    vals = range(Rmin, Rmax, 1000)
    s = Int(sqrt(length(coefs)))
    cx = [range(Rmin, Rmax, s);]
    cy = [range(Rmin, Rmax, s);]
    scale = rbfscale(cx, cy, k)
    mat = [interpolant(rbf, Rfrom, Rto, coefmat(coefs), cx, cy, scale)
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
                     model_kwargs = (; ndc = 4, ngc = 4, normalize = false),
                     optim_kwargs = (; ad = AD));


AD = AutoEnzyme()
AD = ADTypes.AutoForwardDiff()
mdl = norm(data; densscale = 1.5, ndc = 4, ngc = 4, normalize = false);
mdl = dist(data; densscale = 1.5, ndc = 4, ngc = 4, normalize = false);
mdl = distdense(data; densscale = 1.5, ndc = 4, ngc = 4, normalize = false);

mdl = distgeo(data; densscale = 1.5, ndc = 4, ngc = 4, normalize = false);
chn = Turing.sample(mdl.mdl, NUTS(100, .7), 10, progress = true)


plot(chn[:lp], label = "$(round(chn[:lp][end], digits = 2))")
denscoefs = extract_coefs(chn[end, :, 1], "ζ")
geocoefs = extract_coefs(chn[end, :, 1], "η")

preds = returned(mdl.mdl, chn[end])[1];
df = DataFrame(; data.df.fromdist, data.df.todist, data.df.flows, preds, data.df.dist);
plot(plotfit(df.flows, df.preds),
     title = "Cor: $(corround(log.(df.flows), log.(df.preds)))")

plotdist(df, :preds)

net = calc_net_df(df);
plot(plotnet(net), title = "cor: $(corround(net.nmr, net.nmrp))")

districts = year(CSV.read("../data/districts.csv", DataFrame), 2017)
mat = heat(denscoefs, districts);

w = zeros(4, 4)
w[3, 3] = 2
w[1, 1] = 1
w[4, 1] = 1

coefmat(w)
heat(w, 1.0, districts);

plotgeo(w, districts)
