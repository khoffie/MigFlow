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

function fit_map(; model, inits, lower, upper, algo, iters, reltol, dt)
    @printf("Algorithm = %s\n", nameof(typeof(algo)))
    @printf("Number of iterations = %.f\n", iters)
        
    fit = maximum_a_posteriori(model, algo; 
    adtype = AutoReverseDiff(), 
    initial_params = inits, lb = lower, ub = upper,
    maxiters = iters, maxtime = 600, reltol = reltol, 
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

function write_out(; mod_name, opts, preds)
    CSV.write("./fitted_models/opti$mod_name.csv", opts)
    CSV.write("./predictions/FlowDataPreds$mod_name.csv", preds)
    @printf("Predicions and optims saved for model %s\n", mod_name)
end    

function create_model3(; dt, dists, ncoefs, flow_th)
    dists = @orderby(dists, levelcode.(dists.distcode)) ## make sure the district dataframe is sorted by the level code of the dists
    distdens = dists.density
    distdens = distdens ./ maximum(distdens)
    distdens = distdens .- mean(distdens)
    ## ditrict density on a scale definitely between -1 and 1 most likely more like -0.5, 0.5 but not exactly

    dt2 = dt[in.(dt.fromdist, Ref(dists.distcode)) .&& in.(dt.todist, Ref(dists.distcode)), :]
    droplevels!(dt2.fromdist)
    droplevels!(dt2.todist)
    dt2 = dt2[dt2.flows .> flow_th, :]

    Ndist = length(levels(dt2.fromdist))
    Nages = length(levels(dt2.agegroup))
    popgerm = 73000
    meddist = 293.0 
    netactual = calcnet(dt2.flows,
                        levelcode.(dt2.fromdist),
                        levelcode.(dt2.todist),
                        levelcode.(dt2.agegroup),
                        Nages,
                        Ndist)

    model3 = migration3(dt2.flows, sum(dt2.flows), levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                        dt2.frompop, dt2.topop, popgerm, dt2.distance,
                        levelcode.(dt2.agegroup),
                        Nages,
                        dists.xcoord, dists.ycoord, distdens,dists.pop,
                        Ndist, meddist, netactual, ncoefs)
    return (model3, dt2)
end
