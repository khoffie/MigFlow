using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random, ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta, DataFrameMacros,
    StatProfilerHTML
includet("debughelpers.jl")
Random.seed!(20240719)

# munis = CSV.read("./data/munis_pop.csv", DataFrame)
# rename!(munis,Dict(:age_group => :agegroup))
coords_dt = CSV.read("./data/districts.csv", DataFrame)
sims = true
if sims
#    dt = CSV.read("data/simulations.csv",DataFrame)
#    rename!(dt,Dict("predict" => "flows"))
#    dt.flows = round.(Int32,dt.flows)
    dt,knowndesir = simdatafromtemplate(Xoshiro(20240725),"data/flowtemplate.csv")
    rename!(dt,Dict(:age_group => :agegroup,:distance => :dist,:flow => :flows))
else
    dt = CSV.read("/home/donkon/Documents/GermanMigration/data/FlowDataGermans.csv", DataFrame)
end

dt = CSV.read("/home/donkon/Documents/GermanMigration/data/FlowDataGermans.csv", DataFrame)
dt.fromdist = categorical(dt.fromdist)
dt.todist = categorical(dt.todist)
dt.agegroup = categorical(dt.agegroup)
levels!(dt.agegroup,["below18","18-25","25-30","30-50","50-65","above65"])
rename!(dt, Dict(:dist => :distance))

includet("fitmodel3.jl")

x = ["a", "b"]
x_cat = categorical(x)

levelcode(x_cat)
