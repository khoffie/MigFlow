using Turing, SliceSampling, StatsPlots, LinearAlgebra, LogDensityProblems

include("TemperedModel.jl")


@model function foo()
    a ~ MvNormal(fill(0,3),1.0)
end

sam1 = sample(foo(),externalsampler(SliceSampling.HitAndRun(SliceSteppingOut(2.))),MCMCThreads(),10,4,
    initial_params=fill(fill(10.0,3),4))

plot(sam["a[1]"])
plot(sam[:lp])


l = logjoint(foo(),sam)

sam.value[:,4,:] .= l
plot(sam[:lp])


@model function foo2()
    a ~ MixtureModel([MvNormal([0,0],1),MvNormal([4,4],1)],[0.8,0.2])
end
modval = foo2()

#sam = sample(Turing.LogDensityFunction(modval),TemperedSampler(MH(0.1^2*I(2)),[1.0,.8,.6]),1000,
#    initial_params=[1.0,1.0])
#plot(sam)


dtm = TemperedModel(foo(),1.2)

sam = sample(dtm,SliceSampling.HitAndRun(SliceSteppingOut(0.5)),MCMCThreads(),10,4,
    initial_params=fill(fill(10.0,3),4))


ch = Chains([sam[i][j].params[k] for j in eachindex(sam[1]), k in eachindex(sam[1][1].params), i in eachindex(sam)])

parms = [sam[i][j].params[k] for j in eachindex(sam[1]), k in eachindex(sam[1][1].params), i in eachindex(sam)]

lp =  [sam[i][j].lp for j in eachindex(sam[1]),k in 1:1, i in eachindex(sam)]

hcat(parms,lp)


