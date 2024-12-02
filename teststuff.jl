using Revise
includet("main.jl")

path = "./manuscript_input/144temp"
fchains = chainnames(path, "144")
chain = deserialize(joinpath(path, fchains[1]))
plot(chain[:lp])
