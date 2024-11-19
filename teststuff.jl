includet("src/datafitting.jl")

path = "./manuscript_input/tempered"
f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
geodat, agedat = loadallGermData(sample = false, positive_only = true)
agedat = filter(:agegroup => ==("30-50"), agedat)

p = densitychains(chain, agedat, geodat, 1000)
p
p2 = densitychains(chain, agedat, geodat, 1000)
plot(p, p2)
