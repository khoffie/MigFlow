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

gen_random_inits = function(Nages, ncoefs)
    inits = [
        fill(0.0, Nages); #a
        fill(0.25, Nages); #b
        fill(2.0, Nages); #c
        fill(2.0, Nages); #d0
        [.25, 3.0, 5.0]; #dscale, neterr and logconst
        fill(0.0, Nages); #kd
        fill(0.0, Nages * ncoefs) #desire
        ]
  return inits
end

gen_bounds = function(Nages, ncoefs, cheby_lb, cheby_ub)
    a_lb = - 5.5
    a_ub = 5.5
    b_lb = 0.0
    b_ub = 200.0
    c_lb = .05
    c_ub = 10
    d0_lb = 0.0
    d0_ub = 50.0
    ne_lb = .5
    ne_ub = 10.0
    lc_lb = -30.0
    lc_ub = 30.0
    kd_lb = -2.5
    kd_ub = 2.5
        
    lower = [fill(a_lb, Nages); 
            fill(b_lb, Nages); 
            fill(c_lb, Nages); 
            fill(d0_lb, Nages); 
            [0.02]; # dscale in fractions of meddist
            [ne_lb, lc_lb];
            fill(kd_lb, Nages); 
            cheby_lb * ones(ncoefs * Nages)]

    upper = [fill(a_ub, Nages); 
            fill(b_ub, Nages); 
            fill(c_ub, Nages); 
            fill(d0_ub, Nages); 
            [3.0]; #dscale in fractions of meddist
            [ne_ub, lc_ub];
            fill(kd_ub, Nages); 
            cheby_ub * ones(ncoefs * Nages)]
    return(lower, upper)
end

fit_map = function(model, inits, lower, upper, iters, dt)
    fit = maximum_a_posteriori(model, LBFGS(); 
    adtype = AutoReverseDiff(), 
    initial_params = inits, lb = lower, ub = upper,
    maxiters = iters, maxtime = 600, reltol = .05, 
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
