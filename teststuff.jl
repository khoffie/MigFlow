using LogDensityProblems: logdensity
using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise

path = "./manuscript_input/tempered"
f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
flows = CSV.read("./data/FlowDataGermans.csv", DataFrame)
districts = CSV.read("./data/districts.csv", DataFrame)

p = densitychains(chain, flows, 10, districts)
p
x :: Chains = chain
