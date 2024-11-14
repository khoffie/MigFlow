using LogDensityProblems: logdensity
using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise

path = "./manuscript_input/1000slice"


chain = deserialize(joinpath(path, "germchain_2017_30-50.csv"))
geodat, agedat = loadallGermData(sample = false, positive_only = true)

modl = usmodel(agedat.flows, sum(agedat.flows),
    levelcode.(agedat.fromdist), levelcode.(agedat.todist),
    median(geodat.pop), agedat.dist,
    geodat.xcoord, minimum(geodat.xcoord), maximum(geodat.xcoord),
    geodat.ycoord, minimum(geodat.ycoord), maximum(geodat.ycoord),
    geodat.logreldens, minimum(geodat.logreldens), maximum(geodat.logreldens),
    geodat.pop, nrow(geodat), 100.0, 36, 36, true) ## nothing == neta

logjoint(modl, chain)



logprobs = zeros(nsamples)




logdensity(ld, chain.value[i, 1 : end -1, j])
logprobs

chain[: , :lp, :]


?logdensity
weird = -3.2828078659334556

for i in eachindex(modl.args)
    println(sum(weird .== modl.args[i]))
end


path
chain

@model function foo()
    a ~ MvNormal(fill(0,3),1.0)
end

sam = Turing.sample(foo(),
                    externalsampler(SliceSampling.HitAndRun(SliceSteppingOut(2.))),
                    MCMCThreads(), 10, 4,
                    initial_params = fill(fill(10.0,3), 4))

out = logprob(foo(), sam)
logjoint(foo(), sam)


function logprob(model, chain)
    ld = Turing.LogDensityFunction(model)
    nsamples = size(chain)[1]
    nchains = size(chain)[3]
    out = zeros(nsamples, nchains)
    for j in 1 : nchains
        for i in 1 : nsamples
            out[i, j] = logdensity(ld, chain.value[i, 1 : end -1, j])
        end
    end
    return out
end


sam[:, :lp, :] = out


size(sam)
densmin = 0
densmax = 5000

ncoefs = 36
kd_dist = MvNormal(zeros(ncoefs), fill(40.0, ncoefs))
kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
            rand(kd_dist) ./ 10)
Plots.plot(kdfun, seriestpe = :scatter, alpha = .1)
Plots.heatmap(kdfun)
Plots.surface(kdfun, colorbar=true, ticks=false)


path = "manuscript_input/2024-11-14_14-36-20/"
path = "manuscript_input/2024-11-14_14-52-39"
postprocess(path)

