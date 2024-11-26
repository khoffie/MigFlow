using Revise
includet("main.jl")

path = "./manuscript_input/pretempering"
path = "./manuscript_input/30kMH"
jo = "./writeup/juliaout_path.txt"
write(jo, path)

f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
flows, districts = loadallGermData(sample = false, positive_only = true)

postprocess(50, path, false, true)

fromdens, todens, densmin, densmax = densodensd(flows, districts, 10000)
p, v = densheatmap(chain, fromdens, todens, densmin, densmax)
postprocess()

r = 0.1
λ = log(1 - r)
yₜ = 10
y₀ = 10000


decaytime(10000, 1, .1, 10)
