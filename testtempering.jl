using Revise
includet("main.jl")

pos_only = true
districts, flows = loadallGermData(sample = false, positive_only = pos_only)
flows = filter(:agegroup => ==("30-50"), flows)
turingmodel = germmodel(flows, districts, pos_only)

## germdd = (flows = flows, geog = districts, model = turingmodel)
mutable struct germ
    flows :: DataFrame
    districts :: DataFrame
    model :: Any ## better DynamicPPL.Model or TemperedModel
end

germd = germ(flows, districts, turingmodel)
vals = gen_inits()
vals.optis = vals.inits
path = "./manuscript_input/temperedtest/"
chainout = joinpath(path, "germchain_2017_30-50.csv")

sampler = SliceSampling.HitAndRun(SliceSteppingOut(0.25))

germd, vals, chain = runsampling(TemperedModel(germd.model, 1000.0),
                                 germd, sampler, vals, chainout, 1, 10, 1, printvals = false)


allresults=[]
let temp = 10000.0,
##    inits = vals.optis,
    germd = germd,
    vals = vals,
    chain = Chains([1]),
    restartcount = 0

    global allresults

    while temp > 8000.0
            @label restartsample
        try
            tempmodel = TemperedModel(turingmodel, temp)
            println("Sampling starts for temperature $temp")
            germd, vals, chain = runsampling(tempmodel, germd,
                                             SliceSampling.HitAndRun(SliceSteppingOut(0.25)),
                                             vals, chainout, 4, 10, 1, printvals = false)
        catch e 
            println("Error occurred of type $(typeof(e)) potential restart")
##            println(e)
            if typeof(e) != InterruptException && restartcount < 100
                vals.optis = vals.optis .+ rand(Normal(0.0,0.05),length(vals.optis))
                println(vals.optis[1 : 10])
                restartcount += 1
                @goto restartsample
            else
                break
            end
        end
##        plot(chain[50:100,:lp,:],title="Temperature $temp") |> display
        plot(chain[:lp],title="Temperature $temp") |> display
        push!(allresults,(chain=chain,temp=temp))
        temp = temp * 0.9
        vals.optis = vals.optsam
        vals.optis = vals.optis .+ rand(Normal(0.0,0.05),length(vals.optis))
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
