using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise

path = "manuscript_input/2024-11-06_13-37-41"
chain = deserialize(joinpath(path, "germchain_below18.csv"))

plot(chain[:c])
plot(chain[:lp])

