
"""
dists should be a DataFrame with distcode, pop, density, xcoord, ycoord 
"""
function testmod3(dt, optis, dists, flow_th, map_iters, dovi, dosamp)
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

    opinit = gen_random_inits(Nages, ncoefs)
    lower, upper = gen_bounds(Nages, ncoefs, -10.0, 10.0)
    
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
    ## for size in [.1, .2, .4, 1.0, 2.0, 4.0, 10.0]
    for size in [10.0]
        lower, upper = gen_bounds(Nages, ncoefs, - size, size)
        @printf("Starting Optimization for size = %.2f\n", size) 
        @printf("Chosen MAP iterations = %.f\n", map_iters)
        mapfit, opts, preds = fit_map(model3, opinit, lower, upper, map_iters)
        serialize("data/mapfit3_$size.dat", mapfit)
        CSV.write("./data/opti_model3_$size.csv", opts)
        CSV.write("./data/FlowDataPreds3_$size.csv", preds)        
        opinit = opts
    end
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

    return(fit = mapfit, fitdf = opts, dt2 = preds, samps = fit3, visamps = fit4)
end


