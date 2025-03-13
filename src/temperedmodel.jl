
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

function runtempering(data,
                      sampler,
                      params,
                      inits;
                      outpaths,
                      thinning,
                      start_temp,
                      last_temp,
                      decay = .1,
                      n_samples = 100,
                      final_samples,
                      plt_type = "lp")
    ## runtempering is confusing, because vals is updated in each step
    ## and allresults has the same vals for each iteration. But in
    ## fact these are different as can be seen fro, inspecting the
    ## chains.
    allresults = []
    temp = start_temp
    optis = zeros(length(inits))
    restartcount = 0
    chain = nothing
    while temp > last_temp
        @label restartsample
        try
            tempmodel = TemperedModel(data.model, temp)
            println("Sampling starts for temperature $temp")
            addtemp(x) = replace(x, ".csv" => "temp_$(round(temp, digits = 0)).csv")
            temppaths = Dict([k => addtemp(v) for (k, v) in outpaths])
            if temp * (1 - decay) <= last_temp; n_samples = final_samples; end
            data, optis, chain = runsampling(tempmodel, data,
                                             sampler,
                                             params, fill(inits, 4), chainout = temppaths["chain"],
                                             nchains = 4, nsamples = n_samples, thinning = thinning,
                                             paramtype = "best")
            
        catch e ## this catches all errors but it should only catch the domain error
            println("Error occurred of type $(typeof(e)) potential restart")
            println(e)
            if typeof(e) != InterruptException && restartcount < 5
                inits = inits .+ rand(Normal(0.0,0.05), length(inits))
                println("Restarting with perturped inits")
                println(inits[1 : min(6, length(inits))])
                restartcount += 1
                @goto restartsample
            else
                break
            end
        end
        ##        plot(chain[50:100,:lp,:],title="Temperature $temp") |> display
        if plt_type == "lp"
            display(plot(chain[:lp],title="Temperature $temp"))
        elseif plt_type == "gravity"
            display(gravity(chain, 50, 50, "Temperature $temp"))
        end
        push!(allresults,(chain = chain, temp = temp))
        temp = temp * (1 - decay)
        # perturbance to avoid potential type issue
        inits = optis .+ rand(Normal(0.0, 0.05), length(optis)) 
    end
    return allresults
end
