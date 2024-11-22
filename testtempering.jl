using Revise
includet("main.jl")

pos_only = true
districts, flows = loadallGermData(sample = false, positive_only = pos_only)
flows = filter(:agegroup => ==("30-50"), flows)
turingmodel = germmodel(flows, districts, pos_only)

germd = germ(flows, districts, turingmodel)
vals = gen_inits()
vals.optis = vals.inits
path = "./manuscript_input/temperedtest/"
chainout = joinpath(path, "germchain_2017_30-50.csv")

germd, vals, chain = runsampling(TemperedModel(germd.model, 1000.0),
                                 germd, sampler, vals, chainout, 1, 10, 1, printvals = false)

lastchain = runtempering(germd, vals, 1, 8000)

plot(lastchain[Symbol.([:a,:c,:d0,:e,:dscale,:ktopop,"kd[1]","kd[2]","kd[3]","desirecoefs[1]","desirecoefs[2]","desirecoefs[3]",:lp])])


newinit = lastchain.value[100,1:end-1,1]
newinit = newinit .+ rand(Normal(0,.001),length(newinit))

germd, vals, chain = runsampling(germd, SliceSampling.HitAndRun(SliceSteppingOut(0.25)), vals,
                                    chainout, 4, 100, 10,fill(newinit,4); temp = 4500.0, printvals = false)

serialize("manuscript_input/tempered/dan_temper_test.serial",chain)
plot(chain[Symbol.([:a,:c,:d0,:e,:dscale,:ktopop,"kd[1]","kd[2]","kd[3]","desirecoefs[1]","desirecoefs[2]","desirecoefs[3]",:lp])])
