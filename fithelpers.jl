gen_preds = function(mapmodel, optis)
    ## expects values in col2, names in col1
    chain = Chains([optis[: , 2]], optis[: , 1]) 
    preds = generated_quantities(mapmodel, chain)
    return preds
end

get_params = function(mapfit) 
    params = DataFrame(names = names(mapfit.values, 1), values = mapfit.values.array)
    return params   
end

load_flows = function()
    if ENV["USER"] == "konstantin"
        dt = CSV.read("data/FlowDataGermans.csv", DataFrame)
    elseif ENV["USER"] == "donkon" ## USER on main machine in office
        dt = CSV.read("data/FlowDataGermans.csv", DataFrame)
    elseif ENV["USER"] == "dlakelan"
        dt = CSV.read("data/simulations.csv", DataFrame)
        DataFramesMeta.@transform!(dt,:flows = round.(Int32,:predict),:frompop_ger = :frompop, :topop_ger = :topop)
    end
end    

gen_random_inits = function(Nages, ncoefs)
    inits = [
        rand(Normal(0.0, 1.0), Nages); #a
        rand(Gamma(3.0, 1.0 / 2.0), Nages); #b
        rand(Uniform(1.5, 2.5), Nages); #c
        rand(Uniform(0.0, 0.03), Nages); #d0
        [1.5, 8.0]; # neterr and logisticconst
        fill(0.0, Nages); #kd
        fill(0.0, Nages * ncoefs) # desirecoefs
        ]
  return inits
end

gen_bounds = function(Nages, ncoefs, cheby_lb, cheby_ub)
    kd_lb = -1.5
    kd_ub = 1.5
    c_lb = .05
    c_ub = 5
    d0_lb = 0.0
    d0_ub = .03
        
    lower = [fill(-5.5, Nages); 
            fill(0.0, Nages); 
            fill(c_lb, Nages); 
            fill(d0_lb, Nages); 
            [.05, -30.0];
            fill(kd_lb, Nages); 
            cheby_lb * ones(ncoefs * Nages)]

    upper = [fill(20.0, Nages); 
            fill(20.0, Nages); 
            fill(c_ub, Nages); 
            fill(d0_ub, Nages); 
            [3, 30.0];
            fill(kd_ub, Nages); 
            cheby_ub * ones(ncoefs * Nages)]
    return(lower, upper)
end

fit_map = function(model, inits, lower, upper, iters)
    fit = maximum_a_posteriori(model, BBO_adaptive_de_rand_1_bin() ; adtype = AutoReverseDiff(), 
    initial_params = inits, lb = lower, ub = upper,
    maxiters = iters, maxtime = 600, reltol = .08, 
    progress = true, show_trace = true)    

    opts = DataFrame(names=names(fit.values, 1), 
    values = fit.values.array, 
    inits = inits)
    ## display(opts3)
    display(density(opts.values .- opts.inits))

    chain = Chains([opts[: , 2]], opts[: , 1])
    dt[:, "preds"] = generated_quantities(model, chain)[1][1]
    return(fit, opts, dt)
end
