using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random, ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta, DataFrameMacros,
    StatProfilerHTML
includet("debughelpers.jl")
Random.seed!(20240719)

includet("model3.jl")
##optis = CSV.read("./data/opti_d0_allrows.csv", DataFrame)
optis = CSV.read("./data/opti_d0.csv", DataFrame)

dt = CSV.read("/home/konstantin/Documents/GermanMigration/data/FlowDataGermans.csv", DataFrame)
dt.fromdist = categorical(dt.fromdist)
dt.todist = categorical(dt.todist)
dt.agegroup = categorical(dt.agegroup)
levels!(dt.agegroup,["below18","18-25","25-30","30-50","50-65","above65"])
rename!(dt, Dict(:dist => :distance))

meddist = 293.0  # (or so?)

## Create a districts file which has distcode, pop, density, xcoord, ycoord and save it in the data directory
dists = CSV.read("./data/districts.csv",DataFrame)
dists.distcode = categorical(dists.distcode)


"""
dists should be a DataFrame with distcode, pop, density, xcoord, ycoord 
"""
function testmod3(dt,optis,dists,meddist)
    droplevels!(dists.distcode)
    dists = @sort(dists,levelcode.(dists.distcode)) ## make sure the district dataframe is sorted by the level code of the dists
    distdens = dists.density
    distdens = distdens ./ maximum(distdens)
    distdens = distdens .- mean(distdens)
    ## ditrict density on a scale definitely between -1 and 1 most likely more like -0.5, 0.5 but not exactly
    ncoefs = 64
    Nages = 6 ## inits require it, only later we compute it
    opinit = [optis[:, 2]; [1.5];
              fill(0.0, Nages); rand(Normal(0.0, .4), Nages*ncoefs)]
    lower = [fill(0.0,Nages); fill(0.0,Nages); fill(0.0,Nages); fill(0.0,Nages); [.05];
             fill(-.1, Nages); -40 * ones(ncoefs * Nages)]
    upper = [fill(20.0,Nages); fill(10.0,Nages); fill(5.0,Nages); fill(1.0,Nages); [2];
             fill(.1, Nages); 40.0 * ones(ncoefs * Nages)]

    ##    dt2 = dt[dt.fromdist .in dists.distcode .&& dt.todist .in
    ##    dists.distcode,:]

    ## not working for me. First error: "Whitespace not allowed .in",
    ## after fixing "unexpected comma in array expression"
    dt2 = dt[in.(dt.fromdist, Ref(dists.distcode)) .&& in.(dt.todist, Ref(dists.distcode)), :]
    droplevels!(dt2.fromdist)
    droplevels!(dt2.todist)

    # dt2 = dt[dt.flows .> 0, :]
    dt2 = dt2[dt2.flows .> 0, :]
    
    Ndist = length(unique(dt2.fromdist))
    Nages = length(unique(dt2.agegroup))

    netactual = calcnet(dt2.flows,
                        levelcode.(dt2.fromdist),
                        levelcode.(dt2.todist),
                        levelcode.(dt2.agegroup),
                        Nages,
                        Ndist)

    model3 = migration3(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                        dt2.frompop_ger, dt2.topop, dt2.distance,
                        levelcode.(dt2.agegroup),
                        Nages,
                        dists.xcoord, dists.ycoord, distdens,
                        Ndist, meddist, netactual, ncoefs)
    mapfit3 = maximum_a_posteriori(model3, LBFGS() ; adtype = AutoReverseDiff(), 
                                initial_params = opinit, lb = lower, ub = upper,
                                maxiters = 20, maxtime = 60, reltol = .08)

    opts3 = DataFrame(names=names(mapfit3.values, 1), values=mapfit3.values.array, inits = opinit)


    display(density(opts3.values .- opts3.inits))

   ## @show opinit
    model3_chain = Chains([opts3[: , 2]], opts3[: , 1])
    dt2[:, "preds3"] = generated_quantities(model3, model3_chain)[1][1]

    CSV.write("./data/opti_model3.csv", opts3)
    CSV.write("./data/FlowDataPreds3.csv", dt2)

    (fit = mapfit3, dt2 = dt2)
end

# testmod3(dt, optis, dists, meddist)

## try it out:

## smallerdists = @subset(dists,dists.density .< median(dists.density))
smallerdists = dists[dists.density .< 0.5 * median(dists.density), :]
# testmod3(dt,optis,smallerdists,meddist)

# or run the profiler and we see where the time is being spent:

# using StatProfilerHTML

# @profilehtml testmod3(dt, optis, dists, meddist)
testmod3(dt, optis, dists, meddist)

