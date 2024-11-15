
using Pkg
Pkg.activate(".")
using CSV, DataFrames, FixedWidthTables, DataFramesMeta, CategoricalArrays,  LibGit2
using StatsBase, StatsFuns, StatsPlots, Distributions, Random, LaTeXStrings, Plots
using Turing, ReverseDiff, ApproxFun
using Printf, Revise, Dates, Enzyme, Serialization, SliceSampling
using LogDensityProblems, LogDensityProblemsAD, Distances, LinearAlgebra
Enzyme.API.runtimeActivity!(true) ## to deal with an Enzyme bug, per https://discourse.julialang.org/t/enzyme-ready-for-everyday-use-2024/118819/7
import PlotlyJS

includet("models.jl")
includet("samplerows.jl")
#includet("postprocess.jl")
includet("fitandwrite.jl")
includet("TemperedModel.jl")
includet("DataFitting.jl")

function runsampling(alldata, sampler, vals, chainout, nchains, nsamples, thinning, inits; temp, printvals=false)

    println("Sampling starts, temp = $temp")
    ## MH(.1^2*I(length(vals.optis)))
    tem = TemperedModel(alldata.model, temp)
    mhsamp = Turing.sample(tem, sampler, MCMCThreads(),
                           nsamples, nchains, thinning=thinning,
                           initial_params=inits,
                           verbose=true, progress=true)
    mhsamp = make_chains(mhsamp, vals.pars)
    println(mhsamp[:, :lp, 1])
##    mhsamp[:, :lp, :] = logprob(alldata.model, mhsamp)
    Serialization.serialize(chainout, mhsamp)
    println("Sampling finished")
    idx = findmax(mhsamp[:lp][end,])[2]
    vals.optsam = mhsamp.value.data[end, 1:end-1, idx] # last sample, no LP, chain with max LP
    if printvals
        println(vals[[1:10; 43:47], :])
    end
    alldata.flows.preds = generated_quantities(alldata.model, vals.optsam, vals.pars)
    return alldata, vals, mhsamp
end

geodat, agedat = loadallGermData(sample = false, positive_only = true)
agedat.agegroup .== "30-50"
agedat = filter(:agegroup => ==("30-50"), agedat)

vals = gen_inits()
vals.optis = vals.inits

modl = usmodel(agedat.flows, sum(agedat.flows),
    levelcode.(agedat.fromdist), levelcode.(agedat.todist),
    median(geodat.pop), agedat.dist,
    geodat.xcoord, minimum(geodat.xcoord), maximum(geodat.xcoord),
    geodat.ycoord, minimum(geodat.ycoord), maximum(geodat.ycoord),
    geodat.logreldens, minimum(geodat.logreldens), maximum(geodat.logreldens),
    geodat.pop, nrow(geodat), 100.0, 36, 36, true) ## nothing == neta

germd=(flows = agedat, geog = geodat, model = modl)


path = "./manuscript_input/tempered/"
chainout = joinpath(path, "germchain_2017_30-50.csv")

allresults=[]

let temp = 10000.0,
    inits = vals.optis,
    germd = germd,
    vals = vals,
    chain = Chains([1]),
    restartcount = 0

    global allresults
    
    while temp > 4500.0
        
        @label restartsample
        try 
            germd, vals, chain = runsampling(germd, SliceSampling.HitAndRun(SliceSteppingOut(0.25)), vals,
                                 chainout, 4, 100, 1,fill(inits,4); temp = temp, printvals = false)

#postprocess(path, false)
        catch e 
            println("Error occurred of type $(typeof(e)) potential restart")
            println(e)
            if typeof(e) != InterruptException && restartcount < 100
                inits = inits .+ rand(Normal(0.0,0.05),length(inits))
                restartcount += 1
                @goto restartsample
            else
                break
            end
        end

        plot(chain[50:100,:lp,:],title="Temperature $temp") |> display
        push!(allresults,(chain=chain,temp=temp))
        temp = temp * 0.9
        maxlp = findmax(chain[:,:lp,:])
        inits = chain.value[maxlp[2].I[1],1:end-1,maxlp[2].I[2]]
        inits = inits .+ rand(Normal(0.0,0.05),length(inits))
    end
end

lastchain = allresults[end][1]

plot(lastchain[Symbol.([:a,:c,:d0,:e,:dscale,:ktopop,"kd[1]","kd[2]","kd[3]","desirecoefs[1]","desirecoefs[2]","desirecoefs[3]",:lp])])


newinit = lastchain.value[100,1:end-1,1]
newinit = newinit .+ rand(Normal(0,.001),length(newinit))

germd, vals, chain = runsampling(germd, SliceSampling.HitAndRun(SliceSteppingOut(0.25)), vals,
                                    chainout, 4, 100, 10,fill(newinit,4); temp = 4500.0, printvals = false)

serialize("manuscript_input/tempered/dan_temper_test.serial",chain)
plot(chain[Symbol.([:a,:c,:d0,:e,:dscale,:ktopop,"kd[1]","kd[2]","kd[3]","desirecoefs[1]","desirecoefs[2]","desirecoefs[3]",:lp])])
