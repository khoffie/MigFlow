includet("model3.jl")
optis = CSV.read("./data/opti_d0.csv", DataFrame)

ncoefs = 64

age = "18-25"


opinit = [optis[:, 2]; [1.5]; 0 * ones(ncoefs)]
opinit = optis170krows
lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages); [.05]; -4 * ones(ncoefs)]
upper = [fill(20,Nages); fill(10,Nages); fill(5,Nages); fill(1,Nages); [2]; 4 * ones(ncoefs)]

dt2 = dt
## dt2 = dt[dt.agegroup .== age, :]
dt2 = dt[dt.flows .> 1, :]
meddist = median(dt2.distance)
Ndist = length(unique(dt2.fromdist))
Nages = length(unique(dt2.agegroup))


netactual = calcnet(dt2.flows,
                    levelcode.(dt2.fromdist),
                    levelcode.(dt2.todist),
                    levelcode.(dt2.agegroup),
##                    ones(Int64, nrow(dt2)),
                    Nages,
                    Ndist)

model3 = migration3(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop_ger, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup),
                    # ones(Int64, nrow(dt2)),
                    Nages,
                    coords_dt.xcoord, coords_dt.ycoord,
                    Ndist, meddist, netactual, ncoefs)
dt2
mapfit3 = maximum_a_posteriori(model3, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 200, maxtime = 600, reltol = .08)

optis170krows = mapfit3.values

optis_3 = DataFrame(names=names(mapfit3.values, 1), values=mapfit3.values.array)

model3_chain = Chains([optis_3[: , 2]], optis_3[: , 1])
preds3 = generated_quantities(model3, model3_chain)
dt2[:, "preds3"] = generated_quantities(model3, model3_chain)[1][1]

CSV.write("./data/opti_model3.csv", optis_3)
CSV.write("./data/FlowDataPreds3.csv", dt2)

