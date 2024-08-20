using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random,
    ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta,
    StatProfilerHTML, StatsFuns, OptimizationBBO, Printf
includet("debughelpers.jl")
includet("fithelpers.jl")
includet("gen_inits.jl")
includet("models.jl")
includet("fitmodel1.jl")
includet("fitmodel2.jl")
includet("fitmodel3.jl")


Random.seed!(20240719)

dt = load_flows()
## Create a districts file which has distcode, pop, density, xcoord, ycoord and save it in the data directory
dists = CSV.read("./data/districts.csv",DataFrame)
dists.distcode = categorical(dists.distcode)

# testmod3(dt, optis, dists, meddist)

## try it out:

biggerdists = dists[dists.density .> median(dists.density), :]
smallerdists = dists[dists.density .< median(dists.density), :]

sampdists = biggerdists[StatsBase.sample(1:nrow(biggerdists),100; replace=false),:]
sampdists = [sampdists;
            smallerdists[StatsBase.sample(1:nrow(smallerdists),100; replace = false),:]    
            ]


#= choosen_dists = [11000; 14713]            
sampdists = dists[in.(dists.distcode, Ref(choosen_dists)), :]
 =#    


 #= result = testmod3(dt = dt, inits = "", dists = dists,
                        flow_th = -1; map_iters = 1000, 
                        mod_name = "3_alldists_allflows_lbfgs",
                        dosamp = false, dovi = false)
 =#


Nages = 6
ncoefs = 36
flow_th = -1
## not so great
opts_f = "fitted_models/optiNewPopLBFGS.csv"
ib = gen_inits_bounds(Nages = Nages, ncoefs = ncoefs, type = "optimal", 
                      opts_f = opts_f, show =true)

result = @time(testmod3simpl(thedf = dt, dists = dists, 
                             inits = ib[:, "inits"],
                             lowers = ib[:, "lowers"],
                             uppers = ib[:, "uppers"],
                             iters = 1000, preiters = 0, reltol = 1e-2,
                             dosamp = false, dosamptest = false,
                             mod_name = "NewPop", ncoefs = ncoefs, flow_th = flow_th))

#@profilehtml result = testmod3(dt, optis, sampdists, 1, 3, false, false)

# or run the profiler and we see where the time is being spent:

# using StatProfilerHTML

# @profilehtml testmod3(dt, optis, dists, meddist)
#result = testmod3(dt, optis, dists, meddist,false,false)

#chain = Turing.sample(model3, Prior(), 100)
#= 
let 
    thedf = load_flows(),
    dists = dists,
    Ndist = nrow(dists),
    Nages = 6,
    Ncoefs = 36,
    inits = [
        fill(0.0,Nages); #a
        fill(2.0,Nages); #c
        fill(0.1,Nages); #d0
        fill(.2,Nages); #dscale
        [2.0]; #[neterr]
        fill(0.0,Nages); #kd
        fill(0.0,Nages*Ncoefs); #desirecoefs
    ],
    aindx = 1:Nages,
    cindx = aindx .+ Nages,
    d0indx = cindx .+ Nages,
    dscindx = d0indx .+ Nages,
    neterridx = Nages*4+1,
    logistindx = Nages*4+2,
    kdidx = Nages*4+3:(Nages*4+3+Nages),
    desiridx = last(kdidx)+1:length(inits),
    lb = inits - .02,
    ub = inits + .02,
    netactual = calcnet(thedf.flows,
                        levelcode.(thedf.fromdist),
                        levelcode.(thedf.todist),
                        levelcode.(thedf.agegroup),
                        Nages,
                        Ndist),
    model3 = migration3(thedf.flows, sum(thedf.flows), levelcode.(thedf.fromdist), levelcode.(thedf.todist),
                        thedf.frompop, thedf.topop, popgerm, thedf.distance,
                        levelcode.(thedf.agegroup),
                        Nages,
                        dists.xcoord, dists.ycoord, distdens,dists.pop,
                        Ndist, meddist, netactual, ncoefs)



    lb[aindx] .= inits[aindx] .- 2.0
    ub[aindx] .= inits[aindx] .+ 2.0
    lb[logistindx] = inits[logistindx] .- 3.0
    ub[logistindx] = inits[logistindx] .+ 3.0

    ## first make a values be approximately correct
    vals = maximum_a_posteriori(model3,BBO_adaptive_de_rand_1_bin_radiuslimited(),
                        init_params=inits,lb=lb,ub=ub,maxiters = 50, maxtime = 120, progress=true)
    inits = vals.value.array

    # lock those values in place:
    lb[aindx] .= inits[aindx] .- .05
    ub[aindx] .= inits[aindx] .+ .05
    lb[logistindx] = inits[logistindx] .- 0.05
    ub[logistindx] = inits[logistindx] .+ 0.05

    # relax c, dscale, d0, kd
    lb[cindx] .= 0.25 
    ub[cindx] .= 5 
    lb[dscindx] .= .1
    ub[dscindx] .= 2.0
    lb[d0indx] .= .02
    ub[d0indx] .= 3.0
    lb[kdidx] .= -10.0
    ub[kdidx] .= 10.0

    vals = maximum_a_posteriori(model3,BBO_adaptive_de_rand_1_bin_radiuslimited(),
                    init_params=inits,lb=lb,ub=ub,maxiters = 50, maxtime = 120, progress=true)
    inits = vals.value.array

    serialinits = copy(inits)
    serialize("fitted_models/serial_init_finding.dat",serialinits)

    ## narrow the window on c and dscale and d0
    lb[cindx] .= inits[cindx] .- .1
    ub[cindx] .= inits[cindx] .- .1
    lb[dscindx] .= inits[dscindx] .- .05
    ub[dscindx] .= inits[dscindx] .+ .05
    lb[d0indx] .= inits[d0indx] .- .01
    ub[d0indx] .= inits[d0indx] .+ .01
    lb[kdidx] .= inits[kdidx] .- 1.0
    ub[kdidx] .= inits[kdidx] .+ 1.0

    ## winden the box on the chebyshev coefficients
    lb[desiridx] .= -10.0
    ub[desiridx] .= 10.0
    vals = maximum_a_posteriori(model3,BBO_adaptive_de_rand_1_bin_radiuslimited(),
                    init_params=inits,lb=lb,ub=ub,maxiters = 100, maxtime = 120, progress=true)
    inits = vals.value.array
    serialize("fitted_models/serialinits_with_cheby.dat",inits)
    vi(model3,ADVI(10,100; adtype = AutoReverseDiff(true)); θ_init=inits)
end
 =#
