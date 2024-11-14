using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise

path = "manuscript_input/2024-11-06_13-37-41"
chain = deserialize(joinpath(path, "germchain_below18.csv"))

plot(chain[:c])
plot(chain[:lp])

path = readline("./writeup/juliaout_path.txt")

readdir(path)

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
predict(modl, chain[:, :, 1])

path
chain

@model function foo()
    a ~ MvNormal(fill(0,3),1.0)
end

sam = Turing.sample(foo(),
                    externalsampler(SliceSampling.HitAndRun(SliceSteppingOut(2.))),
                    MCMCThreads(), 10, 4,
                    initial_params = fill(fill(10.0,3), 4))

chain


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

