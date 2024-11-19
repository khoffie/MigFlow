using LogDensityProblems: logdensity
using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise

path = "./manuscript_input/slice1000"
f = "germchain_2017_30-50.csv"
chain = deserialize(joinpath(path, f))

p1 = densitychains(chain; densmin = -50, densmax = 50)
p2 = densitychains(chain; densmin = -2.5, densmax = 2.5)

plot(p1, p2)
p2
