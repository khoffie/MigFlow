
struct TemperedModel{M}
    model::M
    temperature::Float64 
end

TemperedModel(m::M,t::Float64) where M <:DynamicPPL.Model = TemperedModel(Turing.LogDensityFunction(m),t) 


function LogDensityProblems.logdensity(m::TemperedModel,params)
    lp = LogDensityProblems.logdensity(m.model,params)
    return lp/m.temperature
end

function LogDensityProblems.capabilities(m::TemperedModel)
    LogDensityProblems.capabilities(m.model)
end

function LogDensityProblems.dimension(m::TemperedModel)
    LogDensityProblems.dimension(m.model)
end

function LogDensityProblems.logdensity_and_gradient(m::TemperedModel,p)
    (l,g) = logdensity_and_gradient(m.model,p)
    return (l/m.temperature,g/m.temperature)
end

function make_chains(slices,parmnames)
    nchains = length(slices) # what if there's only one chain? this might not be a length 1 array
    parms = [slices[i][j].params[k] for j in eachindex(slices[1]), k in eachindex(slices[1][1].params), 
        i in eachindex(slices)]
    lp =  [slices[i][j].lp for j in eachindex(slices[1]),k in 1:1, i in eachindex(slices)]
    ch = Chains(hcat(parms,lp), vcat(parmnames,[:lp]))
end


function testtempered()
    @model function foo()
        a ~ MvNormal([0,0,0], 1)
    end
    tem = TemperedModel(foo(), 1.5)
    nchains = 1
    sam = Turing.sample(tem, SliceSampling.HitAndRun(SliceSteppingOut(0.5)),
                        MCMCThreads(), 100, nchains, initial_params=fill(fill(1.0, 3), nchains))
    chain = make_chains(sam, ["a1","a2","a3"])
    return chain, sam
end

function runtempering(germd, vals, thinning, temp_th, n_samples = 100)
    allresults=[]
    let temp = 10000.0,
        ##    inits = vals.optis,
        germd = germd,
        vals = vals,
        chain = Chains([1]),
        restartcount = 0

        global allresults

        while temp > temp_th
            @label restartsample
            try
                tempmodel = TemperedModel(germd.model, temp)
                println("Sampling starts for temperature $temp")
                germd, vals, chain = runsampling(tempmodel, germd,
                                                 SliceSampling.HitAndRun(SliceSteppingOut(0.25)),
                                                 vals, chainout, 4, n_samples, thinning, printvals = false)
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
    return lastchain = allresults[end][1]
end
