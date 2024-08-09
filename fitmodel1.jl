includet("model1.jl")
optis = CSV.read("./data/opti_d0.csv", DataFrame)

model1 = migration1(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop_ger, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Nages, Ndist, meddist)

opinit = optis[: , 2]
lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages)]
upper = [fill(20,Nages); fill(10,Nages); fill(5,Nages); fill(1,Nages)]

mapfit1 = maximum_a_posteriori(model1, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 20, maxtime = 60, reltol = .08)

mapfit1.optim_result

opti_params = DataFrame(names=names(mapfit1.values, 1), values=mapfit1.values.array)
CSV.write("./data/opti_d0.csv", opti_params)

model1_chain = Chains([optis[: , 2]], optis[: , 1]) 
dt2[:, "preds"] = generated_quantities(model1, model1_chain)[1]

CSV.write("./data/FlowDataPreds.csv", dt2)
CSV.write("/home/konstantin/Diss/inst/extdata/clean/daniel/FlowDataPreds.csv", dt2)
