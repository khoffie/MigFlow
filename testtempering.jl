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

germd, optis, chain = runsampling(turingmodel, germd, s[:sampler],
                                  vals[!, "params"], inits;
                                  chainout = paths["chain"], nchains = s[:nchains],
                                  nsamples = s[:sample_size], thinning = s[:thinning],
                                  paramtype = "best")

## in runtempering I fill inits for chains inside
results = runtempering(germd, vals[!, "params"], vals[!, "inits"];
                       outpaths = paths, thinning = 1, temp_th = 8500, n_samples = 10)

alldata, p, vals = fitandwritefile(germd, settings, paths)


p = "/home/donkon/Documents/GermanMigration/manuscript_input/temperedtest"
chain = deserialize(joinpath(p, "germchain_2017_50-65.csv"))

plot(chain[:lp])
postprocess(path = "./manuscript_input/temperedtest")

