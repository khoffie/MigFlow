includet("model1.jl")
model1 = migration1(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop_ger, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Ndist, meddist)
opinit = [10; 3; 2; .1]
lower = [0; 0; 0; 0]
upper = [30; 10; 5; 1]

mapfit1 = maximum_a_posteriori(model1, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 200, maxtime = 60, reltol = .008)
mapfit1
initvals = mapfit1.values
