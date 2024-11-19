includet("src/datafitting.jl")

path = "./manuscript_input/tempered"
f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
geodat, agedat = loadallGermData(sample = false, positive_only = true)
agedat = filter(:agegroup => ==("30-50"), agedat)


logreldens = log.(geodat.density / median(geodat.density))
minlrd = minimum(logreldens)
maxlrd = maximum(logreldens)



p = densitychains(chain, flows, districts, 1000)
p
