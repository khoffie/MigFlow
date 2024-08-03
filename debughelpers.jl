using DataFramesMeta,DataFrames

function selectdists(munidata, ndist)
    distdata = @by(munidata,:district, :pop = sum(:pop))
    sortdata = sort(distdata,:pop; rev=true)
    ourdists = sortdata[1:10,:district]
    ourdists = [ourdists ; sample(distdata[11:end,:district],ndist-10; replace = false)]
end

function subsetdists(flows,dists)
    ourflows = flows[flows.fromdist .in dists .&& flows.todist .in dists,:]
    droplevels!(ourflows.fromdist)
    droplevels!(ourflows.todist)
    ourflows
end


function tryvi(flows,modelfun,nsamps,gradsamps,nsteps)
    modl = modelfun(flows.flows,flows.fromdist,flows.todist,
        flows.frompop,flows.topop,flows.dist,flows.agegroup,length(levels(flows.agegroup)),
        length(levels(flows.fromdist)))
    vifit = vi(modl,ADVI(gradsamps,nsteps))
    samp = rand(vifit,nsamps)
    (vifit,samp)
end

## nsamps is how many final samples you want,
## gradsamps is the number of samples used to estimate the ELBO gradient
## nsteps is how many optimizer steps to take to optimize the VI



