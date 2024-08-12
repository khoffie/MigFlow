
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

    ncoefs = 36
    meddist = 293.0  # (or so?)
    Nages = 6 ## inits require it, only later we compute it
    popgerm = sum(dists.pop) # total pop of germay in thousands, used in model

    #=     opinit = [optis[:, 2]; [1.5,-3.0];
                     fill(0.0, Nages); rand(Normal(0.0, .4), Nages*ncoefs)]
 =# 
cheby_lb = - .5
cheby_ub = .5
kd_lb = -1.5
kd_ub = 1.5
c_lb = .05
c_ub = 5
d0_lb = 0
d0_ub = .03
    
opinit = [rand(Normal(0.0, 1.0), Nages); #a
                rand(Gamma(3.0, 1.0 / 2.0), Nages); #b
                rand(Uniform(1.5, 2.5), Nages); #c
                rand(Uniform(d0_lb, d0_ub), Nages); #d0
             [1.5, - 4.0]; #neterr and logisticconst
              fill(0.0, Nages); #kd
              fill(0.0, Nages*ncoefs) # desirecoefs
              ]
              
    lower = [fill(-5.5,Nages); fill(0.0,Nages); fill(c_lb,Nages); fill(d0_lb,Nages); [.05, -10.0];
             fill(kd_lb, Nages); cheby_lb * ones(ncoefs * Nages)]
    upper = [fill(20.0,Nages); fill(20.0,Nages); fill(c_ub,Nages); fill(d0_ub,Nages); [3, 0.0];
             fill(kd_ub, Nages); cheby_ub * ones(ncoefs * Nages)]

    dt2 = dt[in.(dt.fromdist, Ref(dists.distcode)) .&& in.(dt.todist, Ref(dists.distcode)), :]
    droplevels!(dt2.fromdist)
    droplevels!(dt2.todist)

    dt2 = dt2[dt2.flows .> flow_th, :]
    
    Ndist = length(unique(dt2.fromdist))
    Nages = length(unique(dt2.agegroup))

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
    mapfit3 = maximum_a_posteriori(model3, BBO_adaptive_de_rand_1_bin() ; adtype = AutoReverseDiff(), 
                                initial_params = opinit, lb = lower, ub = upper,
                                maxiters = map_iters, maxtime = 600, reltol = .08, 
                                progress = true, show_trace = true)

    opts3 = DataFrame(names=names(mapfit3.values, 1), 
                      values=mapfit3.values.array, 
                      inits = opinit)
    display(opts3)
    display(density(opts3.values .- opts3.inits))

    model3_chain = Chains([opts3[: , 2]], opts3[: , 1])
    dt2[:, "preds3"] = generated_quantities(model3, model3_chain)[1][1]

    CSV.write("./data/opti_model3.csv", opts3)
    CSV.write("./data/FlowDataPreds3.csv", dt2)
    @printf("Chosen MAP iterations = %.f\n", map_iters)
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

    (fit = mapfit3, fitdf = opts3, dt2 = dt2, samps = fit3, visamps = fit4)
end


