includet("model3.jl")
optis = CSV.read("./data/opti_d0.csv", DataFrame)

ncoefs = 64

opinit = [optis[:, 2]; [1.5]; 0 * ones(ncoefs)]
lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages); [.05]; -4 * ones(ncoefs)]
upper = [fill(20,Nages); fill(10,Nages); fill(5,Nages); fill(1,Nages); [2]; 4 * ones(ncoefs)]

dt2 = dt

dt2 = dt[dt[: , 6] .> 1, :]
netactual = calcnet(dt2.flows,
                    levelcode.(dt2.fromdist),
                    levelcode.(dt2.todist),
                    levelcode.(dt2.agegroup),Nages,Ndist)

model3 = migration3(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop_ger, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Nages,
                    coords_dt.xcoord, coords_dt.ycoord,
                    Ndist, meddist, netactual, ncoefs)

dt2
mapfit3 = maximum_a_posteriori(model3, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 20, maxtime = 60, reltol = .08)

mapfit3.values
