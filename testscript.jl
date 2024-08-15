using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random,
    ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta,
    StatProfilerHTML, StatsFuns, OptimizationBBO, Printf
includet("debughelpers.jl")
includet("fithelpers.jl")
includet("gen_inits.jl")
includet("models.jl")
includet("fitmodel1.jl")
includet("fitmodel2.jl")
includet("fitmodel3.jl")


Random.seed!(20240719)

model  = "3401dists-1flow_th10itersBBO_adaptive_de_rand_1_bin"

optis = CSV.read("./fitted_models/opti" * model * ".csv", DataFrame)


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

sampdists = biggerdists[StatsBase.sample(1:nrow(biggerdists),100; replace=false),:]
sampdists = [sampdists;
            smallerdists[StatsBase.sample(1:nrow(smallerdists),100; replace = false),:]    
            ]


#= choosen_dists = [11000; 14713]            
sampdists = dists[in.(dists.distcode, Ref(choosen_dists)), :]
 =#    


result = testmod3(dt = dt, inits = "", dists = dists, algo = LBFGS(),
                        flow_th = -1; map_iters = 100, 
                        mod_name = "3",
                        dosamp = false, dovi = false)

#@profilehtml result = testmod3(dt, optis, sampdists, 1, 3, false, false)

# or run the profiler and we see where the time is being spent:

# using StatProfilerHTML

# @profilehtml testmod3(dt, optis, dists, meddist)
#result = testmod3(dt, optis, dists, meddist,false,false)

#chain = Turing.sample(model3, Prior(), 100)

#check_inits(6, 1)

