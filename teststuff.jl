using Revise
includet("main.jl")

path = "./manuscript_input/30kMH"
jo = "./writeup/juliaout_path.txt"
write(jo, path)

f = "germchain_2017_30-50.csv"

path = "./manuscript_input/30kMH"
postprocess(path = path, render_doc = false, denstype = "best")

chain = deserialize(joinpath(path, f))
flows, districts = loadallGermData(sample = false, positive_only = true)

path
postprocess(path = path, render_plots = true, denstype = "best")

