using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings

include("src/utils.jl")
include("src/othermodels.jl")
include("src/estimation.jl")
include("src/loadgermdata.jl")
include("src/diag.jl") ## for calc_net_df
include("src/diagplots.jl")
include("src/fullmodel.jl")

districts = CSV.read("data/districts.csv", DataFrame)
addlrd(districts)
df = CSV.read("data/FlowDataGermans.csv", DataFrame)
df = year(age(pos(df), "30-50"), 2017)
df = sample_flows(df, 0.2)
df = joinlrd(df, districts)

mdat = gen_mdat(df, districts; distscale = 100.0, ndc = 1, ngc = 1)
out = estimate(distonly, mdat)
out2 = estimate(gravity, mdat)

df2 = DataFrame(flows = mdat.flows, our = out.preds, gravity = out2.preds,
               dist = mdat.dist)
plotdist(df2, :gravity, 100)
plotdist(df2, :our, 100)
