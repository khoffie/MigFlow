using Revise
includet("main.jl")

# path = "./manuscript_input/144temp"
# fchains = chainnames(path, "144")
# chain = deserialize(joinpath(path, fchains[1]))
# plot(chain[:lp])

# chains = chainnames("/home/donkon/koho/sampling")

# f = "./writeup/juliaout_path.txt"
# write(f, "manuscript_input/completeMH")

# postprocess(render_plots = false)

germ = loadallGermData(0; sample = false, positive_only = true)
fl = germ.flows[germ.flows.year .== 2017, :]
fl = fl[fl.agegroup .== "30-50", :]

mdl = germmodel(fl, germ.geog, "full", true)
germd = (flows = fl, geog = germ.geog, model = mdl)
outpaths = createpaths("./results/testing", "germ", 2017, "30-50")
vals = gen_inits(mdl)
nchains = 4


# sam = MH(.1^2*I(nrow(vals)))
# ## sam = externalsampler(SliceSampling.HitAndRun(SliceSteppingOut(2.)))

# chain = @time(Turing.sample(mdl, sam, MCMCThreads(),
#                            10, 4, thinning = 1,
#                            initial_params = fill(vals.inits, 4),
#                            verbose = true, progress = true))
# # gravity(chain)

# data, values, chains = runsampling(germd.model, germd,
#                                    sam,
#                                    vals.params,
#                                    fill(vals.inits, nchains);
#                                    chainout = outpaths["chain"],
#                                    nchains = nchains,
#                                    nsamples = 10,
#                                    thinning = 1,
#                                    paramtype = "best",
#                                    printvals = false)

fitandwritefile(germd, settings, outpaths)

# chain_full = deserialize("./manuscript_input/completeMH/germchain_2017_30-50.csv")
# chain = deserialize("./results/fundamentals/germchain_2017_30-50.csv")
