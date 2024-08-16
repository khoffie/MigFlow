function gen_preds(mapmodel, optis)
    ## expects values in col2, names in col1
    chain = Chains([optis[: , 2]], optis[: , 1]) 
    preds = generated_quantities(mapmodel, chain)
    return preds
end

function get_params(mapfit) 
    params = DataFrame(names = names(mapfit.values, 1), values = mapfit.values.array)
    return params   
end

function load_flows()
    if ENV["USER"] == "konstantin"
        dt = CSV.read("data/FlowDataGermans.csv", DataFrame)
    elseif ENV["USER"] == "donkon" ## USER on main machine in office
        dt = CSV.read("data/FlowDataGermans.csv", DataFrame)
    elseif ENV["USER"] == "dlakelan"
        dt = CSV.read("data/simulations.csv", DataFrame)
        DataFramesMeta.@transform!(dt,:flows = round.(Int32,:predict),:frompop_ger = :frompop, :topop_ger = :topop)
    end
end    

function fit_map(; model, inits, lower, upper, algo, iters, dt)
    @printf("Algorithm = %s\n", nameof(typeof(algo)))
    @printf("Number of iterations = %.f\n", iters)
        
    fit = maximum_a_posteriori(model, algo; 
    adtype = AutoReverseDiff(), 
    initial_params = inits, lb = lower, ub = upper,
    maxiters = iters, maxtime = 600, reltol = .05, 
    progress = true, show_trace = true)    
    @printf("LP of fit= %.f\n", fit.lp)
    opts = DataFrame(names=names(fit.values, 1), 
    values = fit.values.array, 
    inits = inits)
    ## display(opts3)
    display(density(opts.values .- opts.inits))

    chain = Chains([opts[: , 2]], opts[: , 1])
    dt[:, "preds"] = generated_quantities(model, chain)[1][1]
    return(fit, opts, dt)
end
