includet("model3.jl")
optis = CSV.read("./data/opti_d0.csv", DataFrame)

meddist = 296.0  # (or so?)

## Create a districts file which has distcode, pop, density, xcoord, ycoord and save it in the data directory
dists = CSV.read("./data/districts.csv",DataFrame)


"""
dists should be a DataFrame with distcode, pop, density, xcoord, ycoord 
"""
function testmod3(dt,optis,dists,meddist)

    dists = @sort(dists,levelcode.(dists.distcode)) ## make sure the district dataframe is sorted by the level code of the dists
    distdens = dists.density
    distdens = distdens / maximum(distdens)
    distdens = distdens - mean(distdens)
    ## ditrict density on a scale definitely between -1 and 1 most likely more like -0.5, 0.5 but not exactly
    ncoefs = 64


    opinit = [optis[:, 2]; [1.5]; 0 * ones(ncoefs * Nages)]
    opinit = [optis[:, 2]; [1.5]; rand(Normal(0.0, .4),Nages*ncoefs)]
    lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages); [.05]; -40 * ones(ncoefs * Nages)]
    upper = [fill(20,Nages); fill(10,Nages); fill(5,Nages); fill(1,Nages); [2]; 40 * ones(ncoefs * Nages)]

    dt2 = dt[dt.fromdist .in dists.distcode .&& dt.todist .in dists.distcode,:]
    droplevels!(dt2.fromdist)
    droplevels!(dt2.todist)

    dt2 = dt[dt.flows .> 0, :]

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


    density(opts3.values .- opts3.inits)

    opinit
    model3_chain = Chains([opts3[: , 2]], opts3[: , 1])
    dt2[:, "preds3"] = generated_quantities(model3, model3_chain)[1][1]

    CSV.write("./data/opti_model3.csv", optis_3)
    CSV.write("./data/FlowDataPreds3.csv", dt2)

    (fit = mapfit3, dt2 = dt2)
end



## try it out:

#=

smallerdists = @subset(dists,dists.density .< median(dists.density))

testmod3(dt,optis,smallerdists,meddist)

=#

