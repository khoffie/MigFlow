using DataFramesMeta,DataFrames


"""
This function takes the results of loading data/munis_pop.csv
    and subsets ndist districts using the top 10 by population, then 
    the rest are random
"""
function selectdists(munidata, ndist)
    distdata = @by(munidata,:district, :pop = sum(:pop))
    sortdata = sort(distdata,:pop; rev=true)
    ourdists = sortdata[1:10,:district]
    ourdists = [ourdists ; sample(distdata[11:end,:district],ndist-10; replace = false)]
end


"""
takes a flows dataset, and a collection of districts, and
subsets the flows data to those that only mention the districts
in the collection
"""
function subsetdists(flows,dists)
    ourflows = flows[flows.fromdist .in dists .&& flows.todist .in dists,:]
    droplevels!(ourflows.fromdist)
    droplevels!(ourflows.todist)
    ourflows
end


"""
takes flows and a model constructor function and constructs
a model then does variational inference on it. nsamps
is how many variational samples you want, gradsamps is how many
samples to take to calculate gradients, and nsteps is how many
optimization steps you want to take. Returns the variational
fit itself, and samples from it
"""
function tryvi(flows,modelfun,nsamps,gradsamps,nsteps)
    modl = modelfun(flows.flows,flows.fromdist,flows.todist,
        flows.frompop,flows.topop,flows.dist,flows.agegroup,length(levels(flows.agegroup)),
        length(levels(flows.fromdist)))
    vifit = vi(modl,ADVI(gradsamps,nsteps))
    samp = rand(vifit,nsamps)
    (vifit,samp)
end




