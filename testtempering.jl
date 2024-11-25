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

results = runtempering(germd, vals[!, "params"], vals[!, "inits"];
                       outpaths = paths, thinning = 1, temp_th = 8500, n_samples = 10)

fill(retparams(results[1].chain, "best"), 4)



results[2].chain[:a]
findmax(results[2].chain[:lp])

results[3].chain[:a]

results
results.vals

alldata, p, vals = fitandwritefile(germd, settings, paths)

s = settings
inits = fill(vals[!, "inits"], s[:nchains])

germd, optis, chain = runsampling(turingmodel, germd, s[:sampler],
                                  vals[!, "params"], inits;
                                  chainout = paths["chain"], nchains = s[:nchains],
                                  nsamples = s[:sample_size], thinning = s[:thinning],
                                  paramtype = "best")
chain[:lp]

function runsampling(model, alldata, sampler, params, inits; chainout, nchains,
                     nsamples, thinning, printvals = false)

optis.data
chain


ch = results.chain
maxlp = findmax(ch[:, :lp, :])
vals.optsam = ch.value[maxlp[2].I[1], 1:end-1, maxlp[2].I[2]] ## best overall sample


moreout(germd, paths, vals)
