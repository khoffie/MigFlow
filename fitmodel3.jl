includet("model3.jl")
optis = CSV.read("./data/opti_d0.csv", DataFrame)

ncoefs = 64


opinit = [optis[:, 2]; [1.5]; 0 * ones(ncoefs * Nages)]
opinit = [optis[:, 2]; [1.5]; rand(Normal(0.0, .4),Nages*ncoefs)]
lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages); [.05]; -40 * ones(ncoefs * Nages)]
upper = [fill(20,Nages); fill(10,Nages); fill(5,Nages); fill(1,Nages); [2]; 40 * ones(ncoefs * Nages)]

dt2 = dt
## dt2 = dt[dt.agegroup .== age, :]
dt2 = dt[dt.flows .> 0, :]
meddist = median(dt2.distance)
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
                    coords_dt.xcoord, coords_dt.ycoord,
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

dt2
