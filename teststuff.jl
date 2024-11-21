using Revise
includet("main.jl")
using KernelDensity

path = "./manuscript_input/30kMH"
f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
districts, flows = loadallGermData(sample = false, positive_only = true)
flows = filter(:agegroup => ==("30-50"), flows)

p = densitychains(chain, flows, districts, 1000, 0.01, "Density preferences, 30-50, 2017")


fromdens, todens, densmin, densmax = densodensd(flows, districts, 10000)
p, v = densheatmap(chain, fromdens, todens, densmin, densmax)


