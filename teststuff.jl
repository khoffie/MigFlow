using Revise
includet("main.jl")

path = "./manuscript_input/30kMH"
jo = "./writeup/juliaout_path.txt"
write(jo, path)

f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
flows, districts = loadallGermData(sample = false, positive_only = true)

postprocess(50, path, false, true)



