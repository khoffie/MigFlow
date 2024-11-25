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

results = runtempering(germd, vals; outpaths = paths, thinning = 1, temp_th = 8000, n_samples = 10)
results.vals

fitandwritefile(germd, settings, paths)

results[end].vals

ch = results.chain
maxlp = findmax(ch[:, :lp, :])
vals.optsam = ch.value[maxlp[2].I[1], 1:end-1, maxlp[2].I[2]] ## best overall sample


moreout(germd, paths, vals)
