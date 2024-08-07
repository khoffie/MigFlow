includet("model1.jl")

### one age group at a time
dt2 = dt[dt.agegroup .== "below18", :]
dt2 = dt
model1 = migration1old(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop_ger, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Nages, Ndist, meddist)
opinit = [10; 3; 2; .1]
lower = [0; 0; 0; 0]
upper = [30; 10; 5; 1]
mapfit1 = maximum_a_posteriori(model1, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 20, maxtime = 60, reltol = .08)
mapfit1
initvals = mapfit1.values


## all age groups together
model1 = migration1new(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop_ger, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Nages, Ndist, meddist)

opinit = [fill(11.0,Nages); fill(3.3,Nages); fill(1.8,Nages); fill(.1,Nages)]
lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages)]
upper = [fill(20,Nages); fill(10,Nages); fill(5,Nages); fill(1,Nages)]

mapfit1 = maximum_a_posteriori(model1, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 200, maxtime = 600, reltol = .08)
mapfit1
initvals = mapfit1.values
