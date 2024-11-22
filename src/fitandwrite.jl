function fitandwritefile(alldata, settings, outpaths)    
    vals = gen_inits()
    vals = runoptim(vals; run = settings[:run_optim], printvals = false)
    alldata, vals = runsampling(alldata, settings[:sampler], vals, outpaths["chain"],
                                settings[:nchains], settings[:sample_size], settings[:thinning];
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
    df = DataFrame(:pars => parnames, :lb => lb, :ub => ub, :inits => ini)
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

# fit = Turing.sample(model3, NUTS(warmup,.8; init_ϵ = 1e-6, 
#                 adtype=AutoReverseDiff(true)), MCMCThreads(), samples, 3,
#                 initial_params = Iterators.repeated(inits), lower = lowers, upper = uppers,    
#                 verbose = true, progress = true)

function runsampling(alldata, sampler, vals, chainout, nchains, nsamples, thinning; printvals=false)
    println("Sampling starts")
    ## MH(.1^2*I(length(vals.optis)))
    mhsamp = Turing.sample(alldata.model, sampler, MCMCThreads(),
                           nsamples, nchains, thinning=thinning,
                           initial_params=fill(vals.optis, nchains),
                           verbose=true, progress=true)
    println("Sampling finished")
    tempered = occursin("TemperedModel", string(alldata.model))
    slice = occursin("HitAndRun", string(sampler))
    if tempered
        println("Make chains from tempered model")
        mhsamp = make_chains(mhsamp, vals.pars)
    end
    if slice && !tempered
        ## lp values are wrong in slice
        println("Compute log probabilities")
        mhsamp[:, :lp, :] = logprob(alldata.model, mhsamp)        
    end
    Serialization.serialize(chainout, mhsamp)
    maxlp = findmax(mhsamp[:, :lp, :])
    vals.optsam = mhsamp.value[maxlp[2].I[1], 1:end-1, maxlp[2].I[2]] ## best overall sample
    if printvals
        println(vals[[1:10; 43:47], :])
    end
    alldata.flows.preds = generated_quantities(alldata.model, vals.optsam, vals.pars)
    return alldata, vals, mhsamp
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
    alldata.geog.desirability = [desirfun(x, y) for (x, y) in zip(alldata.geog.x, alldata.geog.y)]
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
