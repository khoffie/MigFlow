using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random,
    ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta,
    StatProfilerHTML, StatsFuns, OptimizationBBO, Printf,OptimizationNLopt,NLopt

using LogDensityProblems,LogDensityProblemsAD

includet("debughelpers.jl")
includet("fithelpers.jl")
includet("gen_inits.jl")
includet("models.jl")
includet("fitmodel1.jl")
includet("fitmodel2.jl")
includet("fitmodel3.jl")


Random.seed!(20240719)

## Create a districts file which has distcode, pop, density, xcoord, ycoord and save it in the data directory
dists = CSV.read("./data/districts.csv",DataFrame)
dists.distcode = categorical(dists.distcode)

function test(nflow,dists = dists; alg = ParticleSwarm(), niter = 100, nsecs=300, 
                avals = [-1.1,1.1,1.5,-0.6,-1.4,-1.3],
                pctzero = .02)

    Nages = 6
    Ncoefs = 36

    inits = [
        fill(0.0,6); #a
        fill(2.0,6); #c
        fill(1.2,6); #d0
        fill(.75,6); #dscale
        [6.0,150.0]; #[neterr, mm]
        fill(0.0,6); #kd
        fill(0.0,6*Ncoefs); #desirecoefs
    ]

    alldf = load_flows()
    popgerm = sum(dists.pop)
    distdens = dists.density
    distdens = distdens ./ maximum(distdens)
    distdens = distdens .- mean(distdens)
    meddist = median_distance()

    thedf = alldf[StatsBase.sample(1:nrow(alldf),nflow; replace=false),:]
    thedf.rand = rand(Bernoulli(pctzero),nrow(thedf))
    thedf = thedf[thedf.flows .!= 0 .|| thedf.rand .== 1,:]

    @printf("smallest distance: %.2f\n",minimum(thedf.distance))
    @printf("fraction of zeros: %.3f\n",sum(thedf.flows .== 0)/nrow(thedf))
    Ndist = nrow(dists)
    aindx = 1:Nages
    inits[aindx] .= avals

    cindx = aindx .+ Nages
    d0indx = cindx .+ Nages
    dscindx = d0indx .+ Nages
    neterridx = Nages*4+1
    mmindx = neterridx+1
    kdidx = mmindx+1:mmindx+1+Nages
    desiridx = last(kdidx)+1:length(inits)
    lb = inits .- 1.0
    ub = inits .+ 1.0
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
                        Ndist, meddist, netactual, Ncoefs)

    lb[dscindx] .= 0.02
    ub[dscindx] .= 2.0
    lb[d0indx] .= 0.0
    ub[d0indx] .= 10.0
    lb[aindx] .= -4.0
    ub[aindx] .=  4.0
    lb[cindx] .= 1.0
    ub[cindx] .= 3.0
    lb[mmindx] = 1.0
    ub[mmindx] = 500.0
    lb[d0indx] .= 0.0
    ub[d0indx] .= 2.0
    println("""
    a lower bounds are: $(lb[aindx])
    a upper bounds are: $(ub[aindx])
    
    c lower bounds are: $(lb[cindx])
    c upper bounds are: $(ub[cindx])
    """)
    ## first make a values be approximately correct
    println("starting optimize run from\n");
    displayvals(inits)
#    vals = maximum_a_posteriori(model3, alg; adtype = AutoForwardDiff(),
#        initial_params=inits,lb=lb,ub=ub,maxiters = niter, maxtime = nsecs, reltol=1e-5, progress=true)
    vals = maximum_a_posteriori(model3, alg; adtype = AutoForwardDiff(),
                        initial_params=inits,
                        lb=lb,ub=ub,
                        maxiters = niter, maxtime = nsecs, reltol=1e-9, progress=true)
    println("Saving optimum values to fitted_models/serial_init_finding.dat")
    serialize("fitted_models/serial_init_finding.dat",vals)

    println("ended optimize run at\n");
    displayvals(vals.values.array)

    println("Starting variational inference sample")

#    vimod = vi(model3,ADVI(10,100; adtype = AutoReverseDiff(true)))
#    visamp = rand(vimod,100)
    vimod = visamp = nothing
    return vals,model3,vimod,visamp

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


using NamedArrays
function displayvals(vals)
    names = NamedArrays.names(vals)[1] 
    for i in 1:33 
        @show (names[i],vals[i]);
    end
end
function displayvals(vals::Vector{Float64})
    for i in 1:33 
        @show (vals[i]);
    end
end


function runtest()
    #algo = BBO_de_rand_1_bin_radiuslimited()
    #algo = LBFGS()
    #algo = IPNewton();
    #algo = NLopt.LN_NELDERMEAD()
    #algo = NLopt.LN_COBYLA()
    algo = NLopt.LN_BOBYQA()
    init,model,vimod,visamp = test(20000; alg = algo, niter = 500,nsecs = 600,
            pctzero = 1.0); 
    println("plotting fit...."); 
    display(plotfit(init,model))
    displayvals(init.values);
    (init,model,vimod,visamp)
end

function testgradient(model,init)

    vec = init.values.array 
    #vec = vec .+ 0.01 * rand(length(vec))
    ld = LogDensityFunction(model)
    ldg = ADgradient(AutoForwardDiff(),ld)
    (ldval,ldgrad) = LogDensityProblems.logdensity_and_gradient(ldg,vec)
    @show (ldval,ldgrad)

    println("testing gradient at initial condition: ")
    @show init.values.array

    for i in 1:length(vec)
        dx = 1e-11
        dvec = copy(vec)
        dvec[i] = dvec[i] + dx
        lddy = LogDensityProblems.logdensity(ldg,dvec)
        gradi = (lddy - ldval)/dx
        @printf("gradient estimate for dimension %d is %.3f , ADest = %.3f, relerr = %.3f\n",i,gradi,ldgrad[i],(gradi - ldgrad[i])/ldgrad[i])
    end

    println("testing at perturbed condition: ")
    vec = init.values.array .+ rand(Normal(0.0,.001),length(vec))
    (ldval,ldgrad) = LogDensityProblems.logdensity_and_gradient(ldg,vec)
    for i in 1:length(vec)
        dx = 1e-7
        dvec = copy(vec)
        dvec[i] = dvec[i] + dx
        lddy = LogDensityProblems.logdensity(ldg,dvec)
        gradi = (lddy - ldval)/dx
        @printf("gradient estimate for dimension %d is %.3f , ADest = %.3f, relerr = %.3f\n",i,gradi,ldgrad[i],(gradi - ldgrad[i])/ldgrad[i])
    end

end


function testsampling(mod,val)

    Turing.sample(mod,MH(.05^2*I(length(val.values))),MCMCThreads(),1000,4;init_parameters=Iterators.repeated(val.values.array),thinning = 10, discard_initial=300)
#    Turing.sample(mod,NUTS(100,.7),500;init_parameters=val.values.array)
end
