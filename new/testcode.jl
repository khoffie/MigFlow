using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun
using CategoricalArrays, NamedArrays, LaTeXStrings

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
df = year(age(df, "30-50"), 2017)
df = sample_flows(df, 1)
df = joinlrd(df, districts)

mdat = gen_mdat(df, districts; distscale = 100.0, ndc = 1, ngc = 1)
out = estimate(distonly, mdat)
out2 = estimate(gravity, mdat)
out3 = estimate(full, mdat)
out.out

df2 = DataFrame(flows = mdat.flows, our = out.preds, gravity = out2.preds,
               dist = mdat.dist)

df2.res = df2.flows ./ df2.preds

plotdist(df2.flows, df2.our, df2.dist)
plotdist(df2.flows, df2.gravity, df2.dist, true)
