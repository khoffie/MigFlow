using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random,
    ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta,
    StatProfilerHTML, StatsFuns, OptimizationBBO, Printf
includet("debughelpers.jl")
includet("fithelpers.jl")
includet("models.jl")
includet("fitmodel1.jl")
includet("fitmodel2.jl")
includet("fitmodel3.jl")


Random.seed!(20240719)


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

biggerdists = dists[dists.density .> median(dists.density), :]
smallerdists = dists[dists.density .< median(dists.density), :]

sampdists = biggerdists[StatsBase.sample(1:nrow(biggerdists),25;replace=false),:]
sampdists = [sampdists;
            smallerdists[StatsBase.sample(1:nrow(smallerdists),25; replace = false),:]    
            ]

result = @time(testmod3(dt, "",sampdists, 0, 10, false, false))


#= result = testmod3(dt = dt, optis = optis, dists = sampdists, 
                  flow_th = 0, map_iters = 10, false, false)
 =#
#@profilehtml result = testmod3(dt, optis, sampdists, 1, 3, false, false)

# or run the profiler and we see where the time is being spent:

# using StatProfilerHTML

# @profilehtml testmod3(dt, optis, dists, meddist)
#result = testmod3(dt, optis, dists, meddist,false,false)

#chain = Turing.sample(model3, Prior(), 100)
