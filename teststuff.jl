using Revise
includet("main.jl")

path = "./manuscript_input/30kMH"
jo = "./writeup/juliaout_path.txt"
write(jo, path)

f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
flows = CSV.read("./data/FlowDataGermans.csv", DataFrame)
districts = CSV.read("./data/districts.csv", DataFrame)

p = densitychains(chain, flows, 10, districts)
p
x :: Chains = chain

postprocess(50, path, false, true)



