using Revise
includet("main.jl")


vals = gen_inits()
vals.optis = vals.inits
outpath = "./manuscript_input/temperedtest/"
paths = createpaths(outpath, "germ", "2017", "30-50")




pos_only= true
#### works
districts, flows = loadallGermData(sample = false, positive_only = true)
flows = filter(:agegroup => ==("30-50"), flows)
turingmodel = germmodel(flows, districts, true)
germd = (flows = flows, geog = districts, model = turingmodel)
fitandwritefile(germd, settings, paths)
results = runtempering(germd, vals; outpaths = paths, thinning = 1, temp_th = 8000, n_samples = 10)
