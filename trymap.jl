## use optimization to find a good fit
mapfit2 = maximum_a_posteriori(model2, LBFGS() ; adtype = AutoReverseDiff(), 
            initial_params = opinit, maxiters = 20, maxtime = 60, reltol = .08,
            lb = lower, ub = upper)
mapfit2
initvals = mapfit2.values
df = DataFrame(param = names(initvals, 1), estim = values(initvals))
CSV.write("./data/opti_vals.csv", df)
initvals2 = Iterators.Repeated(initvals .+ rand(Uniform(0.0,.50),length(initvals)))
