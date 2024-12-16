using Revise
includet("main.jl")

path = "./manuscript_input/144temp"
fchains = chainnames(path, "144")
chain = deserialize(joinpath(path, fchains[1]))
plot(chain[:lp])

chains = chainnames("/home/donkon/koho/sampling")

f = "./writeup/juliaout_path.txt"
write(f, "./manuscript_input/30kMH")

postprocess(render_plots = false)
