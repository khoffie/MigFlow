"""
dists should be a DataFrame with distcode, pop, density, xcoord, ycoord 
"""
function testmod3(; dt, inits, dists, algo, flow_th, map_iters, mod_name, dovi, dosamp,popgerm,meddist)
    function gen_name(mod_name, Ndists, flow_th, iters, algo, print) 
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
    Ndist = length(levels(dt2.fromdist))
    Nages = length(levels(dt2.agegroup))
    
    # take subsample of flows in long and short distance range, to get the inits started
    longdist = dt2[dt2.distance .> 400.0 .&& rand(Bernoulli(0.1),nrow(dt2)),:]
    shortdist = dt2[dt2.distance .< 400.0 .&& rand(Bernoulli(0.1),nrow(dt2)),:]
    ## better something like (in R) testmod3(..., optis = NULL) and check if optis == NULL
#= 
    if typeof(inits) .== DataFrame
        @printf("Supplied inits used\n")
        inits = inits[:, "values"]
    elseif typeof(inits) .!= DataFrame
        inits = gen_random_inits(Nages, ncoefs)
        @printf("Random inits of length %.f used\n", length(inits))
    end
 =#
    #inits = gen_random_inits(Nages,ncoefs)
    inits = [
        rand(Normal(0.0, 1.0), Nages); #a
        rand(Uniform(2.5, 3.5), Nages); #c
        rand(Uniform(0.0, 0.03), Nages); #d0
        rand(Uniform(0.2, 1), Nages); #dscale
        [1.5, 5.0]; # neterr and logisticconst
        fill(0.0, Nages); #kd
        fill(0.0, Nages * ncoefs) # desirecoefs
        ]

    lower, upper = gen_bounds(Nages, ncoefs)

    model3 = nothing
    netactual = nothing
    mapfit = nothing
    opts = nothing

    for (i,thedf) in zip(1:3,[longdist,shortdist,dt2])
        netactual = calcnet(thedf.flows,
                            levelcode.(thedf.fromdist),
                            levelcode.(thedf.todist),
                            levelcode.(thedf.agegroup),
                            Nages,
                            Ndist)

        model3 = migration3(thedf.flows, sum(thedf.flows), levelcode.(thedf.fromdist), levelcode.(thedf.todist),
                            thedf.frompop, thedf.topop, popgerm, thedf.distance,
                            levelcode.(thedf.agegroup),
                            Nages,
                            dists.xcoord, dists.ycoord, distdens,dists.pop,
                            Ndist, meddist, netactual, ncoefs)

                            ## BBO_adaptive_de_rand_1_bin()


        mod_name = gen_name(mod_name, Ndist, flow_th, map_iters, algo, true)    
        algo2 = algo
        if i == 1
            lower, upper = gen_bounds(Nages, ncoefs, -.02,.02)
            # lock the chebys in place
            lower[1+Nages*4+2:end] .= -.02
            upper[1+Nages*4+2:end] .= .02

            ## precondition the inits by running a BBO optimize
            mapfit, opts, preds = fit_map(model = model3, inits = inits, 
                lower = lower, upper = upper, 
                algo = BBO_adaptive_de_rand_1_bin(), iters = 20, dt = thedf)        
            inits = opts
        elseif i == 2
            lower, upper = gen_bounds(Nages, ncoefs, -.02,.02)
            lower[1+Nages*4+2:end] .= -.02
            upper[1+Nages*4+2:end] .= .02
            ## lock the a values in place
            lower[1:Nages] = inits[1:Nages] - .05*abs.(inits[1:Nages])
            upper[1:Nages] = inits[1:Nages] + .05*abs.(inits[1:Nages])
            #lock logisticconst in place
            lower[1+Nages*4+1] = inits[1+Nages*4+1] - .05*abs.(inits[1+Nages*4+1])
            upper[1+Nages*4+1] = inits[1+Nages*4+1] + .05*abs.(inits[1+Nages*4+1])
            algo = algo2
        else
            lower, upper = gen_bounds(Nages, ncoefs, -1.0,1.0)
        end
        @printf("Number of districts = %.f\n", Ndist)
        @printf("Number of cheby coefs = %.f\n", ncoefs)
        print(show_inits(6, 1))
        print("\n")


        mapfit, opts, preds = fit_map(model = model3, inits = inits, 
                                    lower = lower, upper = upper, 
                                    algo = algo, iters = map_iters, dt = thedf)        
        ##   serialize("data/mapfit3_$size.dat", mapfit)
        CSV.write("./fitted_models/opti$mod_name.csv", opts)
        CSV.write("./predictions/FlowDataPreds$mod_name.csv", preds)        
        inits = opts
        @printf("Model name = %s\n", mod_name)     
    end
    @printf("Begin expanding chebyshev box:\n")
    for size in [.1, .4, 2.0, 4.0, 10.0]
        @printf("Chebyshev box size: %.2f,%.2f\n",-size,size)
        lower, upper = gen_bounds(Nages, ncoefs, -size,size)
        mapfit, opts, preds = fit_map(model = model3, inits = inits, 
                                    lower = lower, upper = upper, 
                                    algo = algo, iters = map_iters, dt = thedf)        
        ##   serialize("data/mapfit3_$size.dat", mapfit)
        mod_name = @sprintf("chebybox%.2f",size)
        CSV.write("./fitted_models/opti$mod_name.csv", opts)
        CSV.write("./predictions/FlowDataPreds$mod_name.csv", preds)        
        inits = opts
        @printf("Model name = %s\n", mod_name)     

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
end


function testmod3simpl(; thedf, dists, iters, preiters)
    dists = @orderby(dists,levelcode.(dists.distcode)) ## make sure the district dataframe is sorted by the level code of the dists
    distdens = dists.density
    distdens = distdens ./ maximum(distdens)
    distdens = distdens .- mean(distdens)
    ## ditrict density on a scale definitely between -1 and 1 most likely more like -0.5, 0.5 but not exactly

    dt2 = dt[in.(dt.fromdist, Ref(dists.distcode)) .&& in.(dt.todist, Ref(dists.distcode)), :]
    droplevels!(dt2.fromdist)
    droplevels!(dt2.todist)
    dt2 = dt2[dt2.flows .> -1, :]

    ncoefs = 36
    Ndist = length(levels(dt2.fromdist))
    Nages = length(levels(dt2.agegroup))
    popgerm = 73000
    meddist = 293.0 
    netactual = calcnet(thedf.flows,
                        levelcode.(thedf.fromdist),
                        levelcode.(thedf.todist),
                        levelcode.(thedf.agegroup),
                        Nages,
                        Ndist)

    model3 = migration3(thedf.flows, sum(thedf.flows), levelcode.(thedf.fromdist), levelcode.(thedf.todist),
                        thedf.frompop, thedf.topop, popgerm, thedf.distance,
                        levelcode.(thedf.agegroup),
                        Nages,
                        dists.xcoord, dists.ycoord, distdens,dists.pop,
                        Ndist, meddist, netactual, ncoefs)

    inits = gen_fixed_inits(Nages, ncoefs)                        
    lowers = gen_bounds(Nages, ncoefs)[1]
    uppers = gen_bounds(Nages, ncoefs)[2]
    @printf("Inits and bounds %s\n", show_inits(Nages, 1))

    mapfit,opts,preds = nothing,nothing,nothing
    if preiters > 0
        ## pre-optimize using a non-gradient optimizer to avoid the worst types of initialization problems
        algo = BBO_adaptive_de_rand_1_bin()
        mapfit, opts, preds = fit_map(model = model3, inits = inits, 
                                        lower = lowers, upper = uppers, 
                                        algo = algo, iters = preiters, dt = thedf)
        inits = opts
        CSV.write("./fitted_models/optiSimpleBBO.csv", opts)
        CSV.write("./predictions/FlowDataPredsSimpleBBO.csv", preds)        
    else
    mapfit, opts, preds = fit_map(model = model3, inits = inits, 
                                    lower = lowers, upper = uppers, 
                                    algo = LBFGS(), iters = iters, dt = thedf)
        CSV.write("./fitted_models/optiSimpleLBFGS.csv", opts)
        CSV.write("./predictions/FlowDataPredsSimpleLBFGS.csv", preds)        
    end        
end
