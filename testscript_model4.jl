using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random,
    ReverseDiff, Revise, RCall, NamedArrays
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta,
    StatProfilerHTML, StatsFuns, OptimizationBBO, Printf,OptimizationNLopt,NLopt
includet("debughelpers.jl")
includet("fithelpers.jl")
includet("gen_inits.jl")
includet("models.jl")
includet("fitmodel1.jl")
includet("fitmodel2.jl")
includet("fitmodel3.jl")


## Random.seed!(20240719)

## Create a districts file which has distcode, pop, density, xcoord, ycoord and save it in the data directory
dists = CSV.read("./data/districts.csv",DataFrame)
dists.distcode = categorical(dists.distcode)

function test(nflow,dists = dists; alg = ParticleSwarm(), niter = 100, nsecs=300, 
                avals = [-1.1,1.1,1.5,-0.6,-1.4,-1.3],
                pctzero = .02, mod_name)

    Nages = 6
    Ncoefs = 36

    inits = [
        fill(0.0, Nages); #a
        fill(2.0, Nages); #c
        fill(1.2, Nages); #d0
        fill(.75, Nages); #dscale
##        fill(600.0, Nages); # mm, per age
    ]

    alldf = load_flows()
    popgerm = sum(dists.pop)
    distdens = dists.density
    distdens = distdens ./ maximum(distdens)
    distdens = distdens .- mean(distdens)
    meddist = median_distance()

    alldf = alldf[alldf.flows .> 0, :]
    thedf = alldf[StatsBase.sample(1:nrow(alldf),nflow; replace=false),:]
    thedf.rand = rand(Bernoulli(pctzero),nrow(thedf))
    thedf = thedf[thedf.flows .!= 0 .|| thedf.rand .== 1,:]
##    thedf = alldf

    @printf("smallest distance: %.2f\n",minimum(thedf.distance))
    @printf("fraction of zeros: %.3f\n",sum(thedf.flows .== 0)/nrow(thedf))
    Ndist = nrow(dists)    

    netactual = calcnet(thedf.flows,
                        levelcode.(thedf.fromdist),
                        levelcode.(thedf.todist),
                        levelcode.(thedf.agegroup),
                        Nages,
                        Ndist)
    model3 = migration4(thedf.flows, sum(thedf.flows), levelcode.(thedf.fromdist), levelcode.(thedf.todist),
                        thedf.frompop, thedf.topop, popgerm, thedf.distance,
                        levelcode.(thedf.agegroup),
                        Nages,
                        dists.xcoord, dists.ycoord, distdens,dists.pop,
                        Ndist, meddist, netactual, Ncoefs)

    aindx = 1:Nages
    cindx = copy(aindx) .+ Nages
    d0indx = copy(cindx) .+ Nages
    dscindx = copy(d0indx) .+ Nages
##    neterrindx = last(dscindx) .+ 1
  ##  mmindx = copy(dscindx) .+ Nages

    lb = inits .- .1
    ub = inits .+ .1
    lb[aindx] .= -4.0
    ub[aindx] .=  4.0
    lb[cindx] .= 1.0
    ub[cindx] .= 3.0
    lb[d0indx] .= 0.0
    ub[d0indx] .= 10.0
    lb[dscindx] .= 0.02
    ub[dscindx] .= 2.0
    # lb[mmindx] .= 500.0
    # ub[mmindx] .= 800.0
    # @show mmindx
    inits[aindx] .= avals
    println("""
    a lower bounds are: $(lb[aindx])
    a upper bounds are: $(ub[aindx])
    
    c lower bounds are: $(lb[cindx])
    c upper bounds are: $(ub[cindx])

    d0 lower bounds are: $(lb[d0indx])
    d0 upper bounds are: $(ub[d0indx])

    dscale lower bounds are: $(lb[dscindx])
    dscale upper bounds are: $(ub[dscindx])

    """)
        # mm lower bounds are: $(lb[mmindx])
    # mm upper bounds are: $(ub[mmindx])

    ## first make a values be approximately correct
    println("starting optimize run from\n");
##     displayvals(inits)
#    vals = maximum_a_posteriori(model3, alg; adtype = AutoForwardDiff(),
    #        initial_params=inits,lb=lb,ub=ub,maxiters = niter, maxtime = nsecs, reltol=1e-5, progress=true)
    println("starting optimization")
    fit = maximum_a_posteriori(model3, alg; adtype = AutoForwardDiff(),
                        initial_params=inits,
                        lb=lb,ub=ub,
                        maxiters = niter, maxtime = nsecs, reltol=1e-9, progress=true)
    opts = DataFrame(names=names(fit.values, 1), 
                    values = fit.values.array, 
                    inits = inits)
    print(opts)
    print("\n")
    chain = Chains([opts[: , 2]], opts[: , 1])
    preds = generated_quantities(model3, chain)
##    print(preds[1:10])
    thedf[:, "preds"] = preds[1]
    write_out(mod_name = mod_name, opts = opts, preds = thedf)

##    serialize("fitted_models/serial_init_finding.dat", )
    return fit, model3

end


function plotfit(vals,model)
    dist = model.args.distance
    flow = Float64.(copy(model.args.flows))
    preds,netflow = generated_quantities(model,vals.values.array,names(vals.values)[1])

    for i in eachindex(flow)
        if flow[i] == 0.0
            flow[i] = preds[i] * rand(LogNormal(-15.0,log(2.0)))
        end
    end
    scatter(dist,log.(flow ./ preds); color = model.args.agegroup,alpha=0.1)
end

function displayvals(vals)
    names = NamedArrays.names(vals)[1] 
    for i in 1:33 
        @show (names[i],vals[i]);
    end
end
function displayvals(vals::Vector{Float64})
    for i in 1:min(length(vals), 33) 
        @show (vals[i]);
    end
end


function runtest()
    algo = BBO_de_rand_1_bin_radiuslimited()
    #algo = LBFGS()
    #algo = IPNewton();
    #algo = NLopt.LN_NELDERMEAD()
    #algo = NLopt.LN_COBYLA()
    algo = NLopt.LN_BOBYQA()
    init,model = test(50000; alg = algo, niter = 5000, nsecs = 6000,
            pctzero = 1.0, mod_name = "esa"); 
    println("plotting fit...."); 
    display(plotfit(init,model))
    savefig("fitted_models/fit.pdf")
    displayvals(init.values);
    (init,model)
end

runtest()
