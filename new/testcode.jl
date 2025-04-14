using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings

include("src/utils.jl")
include("src/othermodels.jl")
include("src/estimation.jl")
include("../src/loadgermdata.jl")
include("src/diag.jl") ## for calc_net_df
include("src/diagplots.jl")
include("src/fullmodel.jl")

districts = CSV.read("../data/districts.csv", DataFrame)
addlrd(districts)
df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
df = year(age(pos(df), "30-50"), 2017)
df = sample_flows(df, 0.2)
df = joinlrd(df, districts)

# plot(plot(out.plt[3], title = "No coefs on frompop and topop"),
#      plot(out2.plt[3], title = "coefs on frompop and topop"),
#      size = (900, 700), layout = (2, 1))

mdat = gen_mdat(df, districts; distscale = 100.0, ndc = 1, ngc = 1)
out = estimate(distonly, mdat)
out2 = estimate(gravity, mdat)

df2 = DataFrame(flows = mdat.flows, our = out.preds, gravity = out2.preds,
               dist = mdat.dist)
plot(plot(plotdist(df2.flows, df2.our, df2.dist), title = "Our model"),
     plot(plotdist(df2.flows, df2.gravity, df2.dist), title = "Gravity"),
     layout = (2, 1))

# dfs = df2[StatsBase.sample(1:nrow(df2), 2000), :]
# plot(xlab = "Distance in km", ylab = L"\log(y / \hat{y})")
# alpha = .2
# hline!([0], color = :darkred, linewidth = 2, label = "")
# scatter!(dfs.dist, log.(dfs.flows ./ dfs.our), label = "Our model",
#         alpha = alpha)
# scatter!(dfs.dist, log.(dfs.flows ./ dfs.gravity), label = "gravity",
#         alpha = alpha)
