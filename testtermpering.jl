function runsampling(alldata, sampler, vals, chainout, nchains, nsamples, thinning; temp, printvals=false)

    println("Sampling starts")
    ## MH(.1^2*I(length(vals.optis)))
    tem = TemperedModel(alldata.model, temp)
    mhsamp = Turing.sample(tem, sampler, MCMCThreads(),
                           nsamples, nchains, thinning=thinning,
                           initial_params=fill(vals.optis, nchains),
                           verbose=true, progress=true)
    mhsamp = make_chains(mhsamp, vals.pars)
    println(mhsamp[:, :lp, 1])
##    mhsamp[:, :lp, :] = logprob(alldata.model, mhsamp)
    Serialization.serialize(chainout, mhsamp)
    println("Sampling finished")
    idx = findmax(mhsamp[:lp][end,])[2]
    vals.optsam = mhsamp.value.data[end, 1:end-1, idx] # last sample, no LP, chain with max LP
    if printvals
        println(vals[[1:10; 43:47], :])
    end
    alldata.flows.preds = generated_quantities(alldata.model, vals.optsam, vals.pars)
    return alldata, vals, mhsamp
end

geodat, agedat = loadallGermData(sample = false, positive_only = true)
agedat.agegroup .== "30-50"
agedat = filter(:agegroup => ==("30-50"), agedat)

vals = gen_inits()
vals.optis = vals.inits

modl = usmodel(agedat.flows, sum(agedat.flows),
    levelcode.(agedat.fromdist), levelcode.(agedat.todist),
    median(geodat.pop), agedat.dist,
    geodat.xcoord, minimum(geodat.xcoord), maximum(geodat.xcoord),
    geodat.ycoord, minimum(geodat.ycoord), maximum(geodat.ycoord),
    geodat.logreldens, minimum(geodat.logreldens), maximum(geodat.logreldens),
    geodat.pop, nrow(geodat), 100.0, 36, 36, true) ## nothing == neta

germd = (flows = agedat, geog = geodat, model = modl)


path = "./manuscript_input/tempered/"
chainout = joinpath(path, "germchain_2017_30-50.csv")

temp = 2^14.

germd, vals, chain = runsampling(germd, SliceSampling.HitAndRun(SliceSteppingOut(1)), vals,
                                 chainout, 4, 100, 1; temp = temp, printvals = false)
postprocess(path, false)
chain[:lp]
