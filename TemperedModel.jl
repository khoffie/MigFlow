
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
        a ~ MvNormal([0,0,0],1)
    end

    f = foo()
    tem = TemperedModel(f,1.5)
    sam = sample(tem,SliceSampling.HitAndRun(SliceSteppingOut(0.5)),MCMCThreads(),10,4,
        initial_params=fill(fill(10.0,3),4))
    return make_chains(sam,["a1","a2","a3"])
end
