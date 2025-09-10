using DataFrames, CSV, CairoMakie, Statistics, StatsBase, Revise
using LaTeXStrings, AlgebraOfGraphics, GeoIO, GeoStats, GLM
using KernelDensity
outp = "/home/konstantin/paper/sections/texdocs/images/"

includet("../src/diagplots.jl")
includet("../src/analyzegeo.jl")
includet("../src/model_helpers.jl")
includet("../src/samplecircle.jl")

df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
df = rmreg(df, :fromdist)
df = rmreg(df, :todist)

di = CSV.read("../data/districts.csv", DataFrame)
di.area = di.pop ./ di.density;

pop = combine(groupby(unique(df, [:fromdist, :year, :agegroup]),
                      [:year, :agegroup]), :frompop => sum => :agepop)
## combine(groupby(pop, :year), :agepop => sum)
shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp")

years = unique(df.year)
ages = unique(df.agegroup)

######################################################################
########################### Distance #################################
res = maincircle(df17, di, false);
dfdist = combine(groupby(df17, :dist),
                 [ nrow => :count,
                   :flows => sum => :flows,
                   [:topop, :frompop] => ((x, y) -> sum(log.(x) + log.(y))) => :pop,
                   [:toarea, :fromarea] => ((x, y) -> sum(log.(x) + log.(y))) => :area
                   ])
res = dropmissing(leftjoin(res, dfdist, on = :radius => :dist))

model(radius, γ, ϵ) = Float64.(radius).^γ .+ ϵ

function plot(res, γ, ϵ)
    res.pred = model(res, γ, ϵ)

    f = Figure(size = (350, 350), fontsize = 10);
    ax1 = Axis(f[1, 1],
              xlabel = "Distance (km)",
              xgridvisible = false, ygridvisible = false)

    xs = StatsBase.sample(res.radius, Weights(res.pop), 10^4);
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :red, label = "Pop")
    xs = StatsBase.sample(res.radius, Weights(res.area), 10^4);
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :blue, label = "Area")
    xs = StatsBase.sample(res.radius, Weights(res.pot), 10^4)
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :purple, label = "a2πr")
    xs = StatsBase.sample(res.radius, Weights(res.count), 10^4)
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :green, label = "Unweighted")
    xs = StatsBase.sample(res.radius, Weights(res.flows), 10^4)
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :yellow, label = "Flows")
    xs = StatsBase.sample(res.radius, Weights(res.count .* res.pred), 10^4)
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :brown, label = "radius^$γ + $ϵ")

    ax2 = Axis(f[1, 2],
               xlabel = "Distance (km)",
               xgridvisible = false, ygridvisible = false)
    res.flows[res.flows .== 0] .= 1

    xs = StatsBase.sample(res.radius, Weights(res.flows ./ res.flows), 10^4);
    lines!(ax2, kde(xs, Gamma(1, 10)), color = :red, label = "Flows / Flows")
    xs = StatsBase.sample(res.radius, Weights(res.flows ./ (res.pred .* res.count)), 10^4);
    lines!(ax2, kde(xs, Gamma(1, 10)), color = :brown, label = "Flows / (radius^$γ + $ϵ)")

    hideydecorations!(ax1)
    hideydecorations!(ax2)
    axislegend(ax1, framevisible = false, patchsize = (5, 10))
    axislegend(ax2, framevisible = false, patchsize = (5, 10))

    titlelayout = GridLayout(f[0, 1], halign = :left, tellwidth = false)
    Label(titlelayout[1, 1], "Movement Distances", halign = :left, fontsize = 15, font = "TeX Gyre Heros Bold Makie")
    Label(titlelayout[2, 1], "2000-2017, all age groups", halign = :left, fontsize = 10)
    rowgap!(titlelayout, 0)

    return res, f
end

r, f = plot(res, (-2.05), 10^(-5))
save(joinpath(outp, "distances.pdf"), f)


df2 = combine(groupby(df, [:fromdist, :todist, :year]),
              :flows => sum => :flows,
              :dist => first => :dist,
              :frompop => sum => :frompop)
df2 = combine(groupby(df2, [:fromdist, :todist]),
              [:flows => mean => :flows,
               :dist => first => :dist,
               :frompop => mean => :frompop])
leftjoin!(df2, year(di, 2017)[!, [:distcode, :pop, :area]], on = :todist => :distcode)

f = Figure(size = (350, 350), fontsize = 10);
ax1 = Axis(f[1, 1],
           xlabel = "Distance (km)",
           xgridvisible = false, ygridvisible = false)

xs = StatsBase.sample(df2.dist, Weights(df2.flows), 10^4)
lines!(ax1, kde(xs, Gamma(1, 10)), color = :red)
xs = StatsBase.sample(df2.dist, Weights(df2.flows_std), 10^4)
lines!(ax1, kde(xs, Gamma(1, 10)), color = :red)


xs = StatsBase.sample(df2.dist, Weights(1 ./ df2.dist.^2), 10^4)
lines!(ax1, kde(xs, Gamma(1, 10)), color = :purple)
xs = StatsBase.sample(df2.dist, Weights(df2.flows ./ (1 ./ df2.dist.^2)), 10^4)
lines!(ax1, kde(xs, Gamma(1, 10)), color = :purple)

xs = StatsBase.sample(df2.dist, Weights(1 ./ df2.dist.^3), 10^4)
lines!(ax1, kde(xs, Gamma(1, 10)), color = :orange)
xs = StatsBase.sample(df2.dist, Weights(df2.flows ./ (1 ./ df2.dist.^3)), 10^4)
lines!(ax1, kde(xs, Gamma(1, 10)), color = :orange)

f

df2.flows_std = df2.flows ./ ((df2.frompop ./ 100) .* (df2.pop ./ 100))
df2
df2.flows_int = Int.(round.(df2.flows, digits = 0))
sort(df2, :flows_std)

mean(df2.flows_std)

@model function mdl(Y, D)
#    α ~ Normal(-5, 1);
    γ ~ Normal(-2, 1)
    ## ϵ ~ Gamma(1, 10)

    λ = γ .* log.(D)
    Y ~ product_distribution(Poisson.(exp.(λ)))
    return λ
end

Turing.sample(mdl(df2.flows_std, df2.dist), NUTS(), 100)
df2
chn2 = Turing.sample(m2, NUTS(), 1000)
res2.preds = returned(m2, chn2[end])[1]
