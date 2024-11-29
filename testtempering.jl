using Revise
includet("main.jl")

outpath = "./manuscript_input/temperedtest/"
paths = createpaths(outpath, "germ", "2017", "30-50")
pos_only= true
#### works
vals = gen_inits()
vals.optis = vals.inits
districts, flows = loadallGermData(sample = false, positive_only = true)
flows = filter(:agegroup => ==("30-50"), flows)
turingmodel = germmodel(flows, districts, true)
germd = (flows = flows, geog = districts, model = turingmodel)

s = settings
inits = fill(vals[!, "inits"], settings[:nchains])

results = runtempering(germd, vals[!, "params"], vals[!, "inits"];
                       outpaths = paths, thinning = 1, temp_th = 8500, n_samples = 10, final_samples = 100)

alldata, p, vals = fitandwritefile(germd, settings, paths)
