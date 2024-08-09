includet("model1.jl")
optis = CSV.read("./data/opti_d0.csv", DataFrame)

Nages = 6
meddist = 293
Ndist = length(unique(dt2.fromdist))

model1 = migration1(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop_ger, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Nages, Ndist, meddist)

opinit = [fill(1,Nages); fill(1,Nages); fill(1,Nages); fill(.1,Nages)]
lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages)]
upper = [fill(20,Nages); fill(10,Nages); fill(5,Nages); fill(1,Nages)]

mapfit1 = maximum_a_posteriori(model1, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 200, maxtime = 600, reltol = .08)

opti_params = DataFrame(names=names(mapfit1.values, 1),
                        values=mapfit1.values.array,
                        old = optis[:, 2])

model1_chain = Chains([optis[: , 2]], optis[: , 1]) 
dt2[:, "preds"] = generated_quantities(model1, model1_chain)[1]

CSV.write("./data/opti_d0_allrows.csv", opti_params[:, 1:2])
CSV.write("./data/FlowDataPreds1.csv", dt2)
