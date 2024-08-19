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

optis = CSV.read("./fitted_models/opti3_alldists_allflows_bbo10iter.csv", DataFrame)

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


 #= result = testmod3(dt = dt, inits = "", dists = dists,
                        flow_th = -1; map_iters = 1000, 
                        mod_name = "3_alldists_allflows_lbfgs",
                        dosamp = false, dovi = false)
 =#


 Nages = 6
 ncoefs = 36
flow_th = -1
## not so great
 opts_f = "fitted_models/wider_kd.csv"
ib = gen_inits_bounds(Nages = Nages, ncoefs = ncoefs, type = "opts", 
                      opts_f = opts_f, show =true)

 result = @time(testmod3simpl(thedf = dt, dists = dists, 
                        inits = ib[:, "inits"],
                        lowers = ib[:, "lowers"],
                        uppers = ib[:, "uppers"],
                        iters = 100, preiters = 0, reltol = 1e-2, dosamp = false, dosamptest = false,
                        mod_name = "NewTest", ncoefs = ncoefs, flow_th = flw_th))

#@profilehtml result = testmod3(dt, optis, sampdists, 1, 3, false, false)

# or run the profiler and we see where the time is being spent:

# using StatProfilerHTML

# @profilehtml testmod3(dt, optis, dists, meddist)
#result = testmod3(dt, optis, dists, meddist,false,false)

#chain = Turing.sample(model3, Prior(), 100)
