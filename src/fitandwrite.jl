function fitandwritefile(alldata, settings, outpaths)    
    vals = gen_inits()
    vals = runoptim(vals; run = settings[:run_optim], printvals = false)
    results = runtempering(alldata, vals, outpaths = outpaths,
                           thinning = 1, temp_th = 8000, n_samples = 10)
    vals.optis = results.vals.optsam
    inits = fill(vals.optis, settings[:nchains])
    alldata, vals.optis, chain = runsampling(alldata.model, alldata, settings[:sampler],
                                             vals.params,
                                             inits, outpaths["chain"], settings[:nchains],
                                             settings[:sample_size], settings[:thinning];
                                             printvals = false)
    moreout(alldata, outpaths, vals)
end

function gen_inits()
    parnames = [["a", "c", "d0", "e", "dscale", "ktopop"];
        ["kd[$i]" for i in 1:36];
        ["desirecoefs[$i]" for i in 1:36]]
    lb = [[-60.0, 0.0, 0.0, 0.0, 1.0, -10.0];
        fill(-50.50, 36);
        fill(-50.50, 36)]
    ub = [[60.0, 20.0, 10.0, 10.0, 15.0, 10.0];
        fill(50.50, 36);
        fill(50.50, 36)]
    ini = rand(Normal(0.0, 0.10), length(ub))
    ini[1:7] .= [-7.6, 1.81, 1.5, 1.5, 5.0, 3.5, 0.0]
    df = DataFrame(:params => parnames, :lb => lb, :ub => ub, :inits => ini)
    return (df)
end

function runoptim(vals; run, printvals=false)
    if run == true
        algo = LBFGS()
        println("Optimization starts")
        mapest = maximum_a_posteriori(alldata.model, algo; adtype=AutoReverseDiff(false),
            initial_params=vals.inits, lb=vals.lb, ub=vals.ub,
            maxiters=50, maxtime=400,
            reltol=1e-3, progress=true)
        println("Optimization finished")
        vals.optis = mapest.values.array
    else
        println("No optimization, using random inits for sampling")
        vals.optis = vals.inits
    end
    if printvals
        println(vals[[1:10; 43:47], :])
    end
    return vals
end

function runsampling(model, alldata, sampler, params, inits; chainout, nchains,
                     nsamples, thinning, printvals = false)
    println("Sampling starts")
    ## MH(.1^2*I(length(vals.optis)))
    chain = Turing.sample(model, sampler, MCMCThreads(),
                           nsamples, nchains, thinning = thinning,
                           initial_params = inits,
                           verbose = true, progress = true)
    println("Sampling finished")
    tempered = occursin("TemperedModel", string(model))
    slice = occursin("HitAndRun", string(sampler))
    if tempered
        println("Make chains from tempered model")
        chain = make_chains(chain, vals.params)
    end
    if slice && !tempered
        ## lp values are wrong in slice
        println("Compute log probabilities")
        chain[:, :lp, :] = logprob(model, chain)        
    end
    Serialization.serialize(chainout, chain)
    maxlp = findmax(chain[:, :lp, :])
    optis = chain.value[maxlp[2].I[1], 1:end-1, maxlp[2].I[2]] ## best overall sample
    # if printvals # rewrite to print DataFrame(inits, optis) or so
    #     println(vals[[1:10; 43:47], :])
    # end
    alldata.flows.preds = generated_quantities(alldata.model, optis, params)
    return alldata, optis.data, chain
end

function moreout(alldata, outpaths, vals)
    (densmin, densmax) = (alldata.model.args.densmin, alldata.model.args.densmax)
    (xmin, xmax) = (alldata.model.args.xmin, alldata.model.args.xmax)
    (ymin, ymax) = (alldata.model.args.ymin, alldata.model.args.ymax)
    kdindx = (1:36) .+ 6 #number of non kd or cheby inits/ priors
    desindx = (1:36) .+ (6 + 36)
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
        vals.optsam[kdindx] ./ 10)
    desirfun = Fun(ApproxFun.Chebyshev(xmin .. xmax) * ApproxFun.Chebyshev(ymin .. ymax),
        vals.optsam[desindx] ./ 10)
    alldata.geog.desirability = [desirfun(x, y) for (x, y) in zip(alldata.geog.xcoord, alldata.geog.ycoord)]
    densvals = range(minimum(alldata.geog.logreldens), maximum(alldata.geog.logreldens), 100)
    densfundf = DataFrame((fromdens=fd, todens=td, funval=kdfun(fd, td)) for fd in densvals, td in densvals)
    CSV.write(outpaths["geog"], alldata.geog)
    CSV.write(outpaths["densfun"], densfundf)
    CSV.write(outpaths["params"], DataFrame(paramval=vals.optsam, parname=vals.pars))
    CSV.write(outpaths["flows"], alldata.flows)
end

function logprob(model, chain)
    ld = Turing.LogDensityFunction(model)
    nsamples = size(chain)[1]
    nchains = size(chain)[3]
    out = zeros(nsamples, nchains)
    for j in 1:nchains
        for i in 1:nsamples
            out[i, j] = LogDensityProblems.logdensity(ld, chain.value[i, 1:end-1, j])
        end
    end
    return out
end
