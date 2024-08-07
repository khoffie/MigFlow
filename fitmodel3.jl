includet("model3.jl")
coords = CSV.read("./data/district_coords.csv", DataFrame)
levelcode.(dt2.fromdist)
opinit = [mapfit1.values.array; [1]; fill(.1,64)]
lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages); [0.05]; fill(-4,64)]
upper = [fill(20,Nages); fill(10,Nages); fill(5,Nages); fill(1,Nages); [2]; fill(4,64)]

model3 = migration3(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop_ger, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Nages,
                    coords.xcoord, coords.ycoord,
                    Ndist, meddist, netactual, 64)

mapfit3 = maximum_a_posteriori(model3, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 20, maxtime = 60, reltol = .08)

