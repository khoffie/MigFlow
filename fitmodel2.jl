includet("model2.jl")
netactual = calcnet(dt2.flows,
                    levelcode.(dt2.fromdist),
                    levelcode.(dt2.todist),
                    levelcode.(dt2.agegroup),Nages,Ndist)

model2 = migration2(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Nages, Ndist, meddist, netactual)


model2
opinit = [fill(11.0,Nages); fill(3.3,Nages); fill(1.8,Nages); fill(.1,Nages); [1]; 1 .* 100 .* ones(Ndist*Nages)]

opinit = [  [.5, 4.5, 5.7, .53, .2, .14];
            [5.6, 7.1, 7.9, 5.8, 3.8, 3.7];
            [4.0, 4.0, 4.0, 4.5, 4.0, 4.0]; 
            [.6, .8, .75, .64, .5, .47]; 
            [1.0]; 
            1 .* 100 .* ones(Ndist*Nages)]

lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages); [0]; .7 .* 100 .* ones(Ndist*Nages)]
upper = [fill(20,Nages); fill(20,Nages); fill(5,Nages); fill(1,Nages); [2]; 1.3 .* 100 .* ones(Ndist*Nages)]

# appx = vi(model2, ADVI(10, 1000))
# rand(appx, 1000)
## use optimization to find a good fit
##lb = lower, ub = upper,
mapfit2 = maximum_a_posteriori(model2, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 200, maxtime = 60, reltol = .008)


## started 14:59
mapfit2
initvals = mapfit2.values
initvals2 = Iterators.Repeated(initvals .+ rand(Uniform(0.0,.50),length(initvals)))
# df = DataFrame(param = names(initvals, 1), estim = values(initvals))
# CSV.write("./data/opti_vals.csv", df)

fit2 = sample(model2, NUTS(500,.8; adtype=AutoReverseDiff(true)), 100,
              init_params = opinit,
              verbose = true, progress = true)
fit2
plot(fit2)
display(fit2)

mainparms2 = fit2[:,[:a,:b,:c,:d0,:neterr],:]
plot(mainparms2) |> display()


opinit

chain = sample(model2, Prior(), 30)
chain
plot(chain)
loglikelihood(model2, opinit)
