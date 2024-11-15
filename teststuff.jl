using LogDensityProblems: logdensity
using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise


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

