
"""
dists should be a DataFrame with distcode, pop, density, xcoord, ycoord 
"""
function testmod3(; dt, inits, dists, algo, flow_th, map_iters, mod_name, dovi, dosamp)
    gen_name = function(mod_name, Ndists, flow_th, iters, algo, print) 
        algo = nameof(typeof(algo))
        mod_name = string(mod_name) * string(Ndists) * "dists" * string(flow_th) *
                  "flow_th" * string(iters) * "iters" * string(algo)
        if print == true
            @printf("Model name = %s\n", mod_name) 
        end                 
        return mod_name
    end

    droplevels!(dists.distcode)
    dists = @orderby(dists,levelcode.(dists.distcode)) ## make sure the district dataframe is sorted by the level code of the dists
    distdens = dists.density
    distdens = distdens ./ maximum(distdens)
    distdens = distdens .- mean(distdens)
    ## ditrict density on a scale definitely between -1 and 1 most likely more like -0.5, 0.5 but not exactly

    dt2 = dt[in.(dt.fromdist, Ref(dists.distcode)) .&& in.(dt.todist, Ref(dists.distcode)), :]
    droplevels!(dt2.fromdist)
    droplevels!(dt2.todist)
    dt2 = dt2[dt2.flows .> flow_th, :]

    ncoefs = 36
    meddist = 293.0
    Ndist = length(unique(dt2.fromdist))
    Nages = length(unique(dt2.agegroup))
    popgerm = sum(dists.pop) # total pop of germay in thousands, used in model

    ## better something like (in R) testmod3(..., optis = NULL) and check if optis == NULL
    if typeof(inits) .== DataFrame
        @printf("Supplied inits used\n")
        inits = inits[:, "values"]
    elseif typeof(inits) .!= DataFrame
        inits = gen_random_inits(Nages, ncoefs)
        @printf("Random inits of length %.f used\n", length(inits))
    end
    lower, upper = gen_bounds(Nages, ncoefs)
    
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

                        ## BBO_adaptive_de_rand_1_bin()

    mapfit3 = nothing
    opts3 = nothing
    mod_name = gen_name(mod_name, Ndist, flow_th, map_iters, algo, true)    
    ## for size in [.1, .2, .4, 1.0, 2.0, 4.0, 10.0]
    lower, upper = gen_bounds(Nages, ncoefs)
    @printf("Number of districts = %.f\n", Ndist)
    @printf("Number of cheby coefs = %.f\n", ncoefs)
    print(show_inits(6, 1))
    print("\n")
    mapfit, opts, preds = fit_map(model = model3, inits = inits, 
                                lower = lower, upper = upper, 
                                algo = algo, iters = map_iters, dt = dt2)        
    ##   serialize("data/mapfit3_$size.dat", mapfit)
    CSV.write("./fitted_models/opti$mod_name.csv", opts)
    CSV.write("./predictions/FlowDataPreds$mod_name.csv", preds)        
    inits = opts
    
    fit3 = nothing
    
    if dosamp
        fit3 = Turing.sample(model3, NUTS(500,.8; adtype=AutoReverseDiff(true)), 100,
                    init_params = opinit,
                    verbose = true, progress = true)
    end
    fit4 = nothing
    if dovi
        fit4 = Turing.vi(model3,ADVI())
    end
end


