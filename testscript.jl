using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random,
    ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta,
    StatProfilerHTML, StatsFuns, OptimizationBBO
includet("debughelpers.jl")
includet("fithelpers.jl")
includet("models.jl")
includet("fitmodel1.jl")
includet("fitmodel2.jl")
includet("fitmodel3.jl")


Random.seed!(20240719)


##optis = CSV.read("./data/opti_d0_allrows.csv", DataFrame)
#optis = CSV.read("./data/opts1greater0.csv", DataFrame)
optis = CSV.read("./data/opti_d0.csv",DataFrame)


dt = load_flows()
dt.fromdist = categorical(dt.fromdist)
dt.todist = categorical(dt.todist)
dt.agegroup = categorical(dt.agegroup)
levels!(dt.agegroup,["below18","18-25","25-30","30-50","50-65","above65"])
rename!(dt, Dict(:dist => :distance))

## Create a districts file which has distcode, pop, density, xcoord, ycoord and save it in the data directory
dists = CSV.read("./data/districts.csv",DataFrame)
dists.distcode = categorical(dists.distcode)


# testmod3(dt, optis, dists, meddist)

## try it out:

## smallerdists = @subset(dists,dists.density .< median(dists.density))
smallerdists = dists[dists.density .< 0.5 * median(dists.density), :]
smallerdists = dists[dists.density .< 0.5 * median(dists.density), :]

smallerdists = dists[shuffle(1 : nrow(dists))[1:50] , : ]
# testmod3(dt,optis,smallerdists,meddist)

result = testmod3(dt, optis, smallerdists, 1, 100, false, false)


# or run the profiler and we see where the time is being spent:

# using StatProfilerHTML

# @profilehtml testmod3(dt, optis, dists, meddist)
#result = testmod3(dt, optis, dists, meddist,false,false)

#chain = Turing.sample(model3, Prior(), 100)

