using Pkg
Pkg.activate(".")
using CSV, DataFrames, FixedWidthTables, DataFramesMeta, CategoricalArrays
using StatsBase, StatsFuns, StatsPlots, Distributions, Random, StatProfilerHTML
using Turing, OptimizationOptimJL, ApproxFun, OptimizationBBO, OptimizationNLopt, NLopt, ReverseDiff
using Printf, Revise, Dates, Enzyme, Serialization
using LogDensityProblems, LogDensityProblemsAD, Distances, LinearAlgebra
Enzyme.API.runtimeActivity!(true) ## to deal with an Enzyme bug, per https://discourse.julialang.org/t/enzyme-ready-for-everyday-use-2024/118819/7
import PlotlyJS

includet("models.jl")
includet("samplerows.jl")

#using DuckDB


#=
## Try fitting migration models to US data. First we need some US data:
## County migration data is available at:

## https://www.census.gov/data/tables/2020/demo/geographic-mobility/county-to-county-migration-2016-2020.html

## County and state geography (Centroids and areas) are available at:

## https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.2020.html#list-tab-264479560

## County population measures are available at: 

## https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html

=#

#const non48fips = (2,11,15,72) #fips for AK, HI, DC, PR

function getUSflows()
    wd = pwd()
    cd("data")
    try 
        download("https://www2.census.gov/programs-surveys/demo/tables/geographic-mobility/2020/county-to-county-migration-2016-2020/county-to-county-migration-flows/CtyxCty_US.txt","CtyxCty_US_2016-2020.txt")
    finally 
        cd(wd)
    end
end

function countycode(stfips,cofips)
    @sprintf("%d_%d",stfips,cofips)
end


function loadUS48flows(geog = loadUS48geog())
    non48fips = (2,11,15,72) #fips for AK, HI, DC, PR
    flows = FixedWidthTables.read("data/CtyxCty_US_2016-2020.txt", (
        STATEFP = (1:3, String),
        COUNTYFP = (4:6, String),
        RES1YRSTATEFP = (7:9, String),
        RES1YRCOUNTYFP = (10:12,String),
        MOVERSWINFLOW = (382:388,Int32)
        )
        ) |> DataFrame
    flows.RES1YRSTATEFP = tryparse.(Int32,flows.RES1YRSTATEFP)
    flows = flows[.!isnothing.(flows.RES1YRSTATEFP),:]
    flows.STATEFP = tryparse.(Int32,flows.STATEFP)
    flows.COUNTYFP = tryparse.(Int32,flows.COUNTYFP)
    flows.RES1YRCOUNTYFP = tryparse.(Int32,flows.RES1YRCOUNTYFP)
    flows = @chain flows begin
        @subset( .! in.(:STATEFP,Ref(non48fips)) .&& .! in.(:RES1YRSTATEFP,Ref(non48fips)))
        DataFramesMeta.@transform(@byrow begin :fromcounty = countycode(:RES1YRSTATEFP,:RES1YRCOUNTYFP)
            :tocounty = countycode(:STATEFP,:COUNTYFP)
        end )
    end
    flows.fromcounty = categorical(flows.fromcounty)
    flows.tocounty = categorical(flows.tocounty)
    levels!(flows.fromcounty,geog.countyid)
    levels!(flows.tocounty,geog.countyid)
    rename!(flows,:MOVERSWINFLOW => :COUNT)
    select!(flows,[:COUNT,:fromcounty,:tocounty])
    flows
end

function getUSgeog()
    wd = pwd()
    cd("data")
    try 
        download("https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2020_Gazetteer/2020_Gaz_counties_national.zip","2020_Gaz_counties_national.zip")
        download("https://www2.census.gov/geo/docs/reference/codes2020/national_county2020.txt","national_county2020.txt") # county identifiers
        run(`unzip 2020_Gaz_counties_national.zip`)
        out = open("2020_Gaz_counties_national.tsv","w")
        for l in eachline("2020_Gaz_counties_national.txt")
            l = strip(l) # get rid of spurious spaces
            println(out,l)
        end
    finally 
        cd(wd)
    end
end

function loadUSgeog()
    geog = CSV.read("data/2020_Gaz_counties_national.tsv",DataFrame; delim="\t")
    countyid = @select(CSV.read("data/national_county2020.txt",DataFrame),:STATE,:STATEFP,:COUNTYFP,:COUNTYNS)
    geog = leftjoin(geog,countyid,on = :ANSICODE => :COUNTYNS)
    geog.countyid = categorical(countycode.(geog.STATEFP,geog.COUNTYFP))
    sort!(geog,[:countyid])
end

function loadUS48geog()
    g = loadUSgeog()
    g = g[.! in.(g.STATE,Ref(("AK","HI","PR","DC"))),: ]
    droplevels!(g.countyid)
    sort!(g,[:countyid])
    g
end


function getUScountypop()
    wd = pwd()
    cd("data")
    try
        download("https://www2.census.gov/programs-surveys/popest/datasets/2010-2020/counties/totals/co-est2020-alldata.csv","co-est2020-alldata.csv")
    finally
        cd(wd)
    end
end



function loadUScountypop()
    countypop = @subset(CSV.read("data/co-est2020-alldata.csv",DataFrame),:COUNTY .!= 0) # filter out state level estimates
    rename!(countypop,Dict(:STATE => :STATEFP,:COUNTY => :COUNTYFP)) # rename the FIPS code columns to indicate they are FIPS and match the geog etc columns
end

function loadGermGeog()
    geog = CSV.read("data/districts.csv",DataFrame)
    geog.distcode = categorical(geog.distcode)
    sort!(geog,[:distcode])
    meddens = median(geog.density)
    geog.logreldens = log.(geog.density ./ meddens)

    geog
end

function loadGermFlows(distcodes)
    fl = CSV.read("data/FlowDataGermans.csv", DataFrame)
    fl = filter(row -> row.fromdist in distcodes, fl)
    fl = filter(row -> row.todist in distcodes, fl)
    fl.fromdist = categorical(fl.fromdist, levels = distcodes)
    fl.todist = categorical(fl.todist, levels = distcodes)
    return fl
end


function loadallGermData(nzeros = 0; sample)
    geog = loadGermGeog()
    if sample == true
        geog = sample_germ(geog)
    end
    flows = loadGermFlows(levels(geog.distcode))
    return (geog=geog, flows=flows)
end

function loadallUSdata(nzeros = 0; sample)

    flows = loadUS48flows()
    geog = loadUS48geog()
    if sample == true
        geog, flows = sample_us(geog, flows)
    end
    
    pop = loadUScountypop()
    pop = DataFramesMeta.@select!(pop,:STATEFP,:COUNTYFP,:POPESTIMATE2016)
    geog = leftjoin(geog,pop,on = [:STATEFP,:COUNTYFP])
    meddens = median(geog.POPESTIMATE2016 ./ geog.ALAND)
    geog.logreldens = log.(geog.POPESTIMATE2016 ./ geog.ALAND / meddens)
    geog.countylevel = levelcode.(geog.countyid)
    sort!(geog,[:countylevel])
    levels!(flows.fromcounty,levels(geog.countyid))
    levels!(flows.tocounty,levels(geog.countyid))

    flowset = Set((f.fromcounty,f.tocounty) for f in eachrow(flows))
    zerosamp = DataFrame((COUNT = 0,fromcounty = f,tocounty = t) for (f,t) in 
        zip(StatsBase.sample(geog.countyid,nzeros),StatsBase.sample(geog.countyid,nzeros)) if !((f,t)  in flowset))
    flows = [flows; zerosamp]

    flows.dist = [haversine((geog[levelcode(r.fromcounty),:INTPTLONG],geog[levelcode(r.fromcounty),:INTPTLAT]), 
                        (geog[levelcode(r.tocounty),:INTPTLONG],geog[levelcode(r.tocounty),:INTPTLAT])) / 1000.0 for r in eachrow(flows)] ## haversine is in m, calculate km
    #netactual::Vector{Float64} = usnetmig(levelcode.(flows.fromcounty),levelcode.(flows.tocounty),flows.COUNT)

    mod = usmodel(flows.COUNT,sum(flows.COUNT),levelcode.(flows.fromcounty),
                levelcode.(flows.tocounty),median(geog.POPESTIMATE2016),flows.dist,
                geog.INTPTLONG,minimum(geog.INTPTLONG),maximum(geog.INTPTLONG),
                geog.INTPTLAT,minimum(geog.INTPTLAT),maximum(geog.INTPTLAT),
                geog.logreldens,minimum(geog.logreldens),maximum(geog.logreldens),
                geog.POPESTIMATE2016,nrow(geog),100.0,36,36)
    (geog=geog, flows = flows, model=mod) #zerosamp=zerosamp,model=mod)
end

function usnetmig(from,to,count)
    nets = zeros(typeof(count[1]),length(unique([from;to])))
    for (f,t,c) in zip(from,to,count)
        nets[f] -= c
        nets[t] += c
    end
    nets
end


function totalflow(fromc, toc, flow, counties)
    fromcindx = levelcode.(fromc)
    tocindx = levelcode.(toc)
    totflow = zeros(length(counties))
    for (fr,to,fl) in zip(fromcindx,tocindx,flow)
        totflow[fr] += fl
        totflow[to] += fl
    end
    totflow
end




function main_interact()
#    densdict = Dict(alldata.county.countyid .=> alldata.county.logreldens)
#    histogram2d([densdict[f] for f in alldata.flows.fromcounty],[densdict[f] for f in alldata.flows.tocounty])

@eval begin
    alldata = loadallUSdata(0) # add no zeros

    #algo = BBO_de_rand_1_bin_radiuslimited()
    algo = LBFGS()
    #algo = NLopt.LN_BOBYQA()

    parnames = [["a","c","d0","dscale","neterr","ktopop"];
        ["kd[$i]" for i in 1:36];
        ["desirecoefs[$i]" for i in 1:36]]


    lb = [[-60.0,0.0,0.0,1.0,0.01,-10.0];
            fill(-50.50,36);
            fill(-50.50,36)]
    ub = [[60.0, 20.0,10.0,15.0,100.0,10.0];
            fill(50.50,36);
            fill(50.50,36)]
    ini = rand(Normal(0.0,0.10),length(ub))
    ini[1:7] .= [-7.6,1.81,1.5,5.0,1.0,3.5,0.0] 
    mapest = Turing.Optimisation.maximum_a_posteriori(alldata.model, algo;
                                                      adtype = AutoReverseDiff(false),
                                                      initial_params = ini,
                                                      lb = lb,ub = ub,
                                                      maxiters = 500, maxtime = 400,
                                                      reltol = 1e-3, progress = true)
    serialize("fitted_models/USmodel_map_$(now()).dat",mapest)
    paramvec = mapest.values.array
    usdiagplots(alldata,paramvec)
    mhsamp = Turing.sample(alldata.model,MH(.1^2*I(length(mapest.values))),100; thinning=50, initial_params = paramvec)
#    mhsamp = Turing.sample(alldata.model,HMCDA(200,.7,1.0; adtype=AutoReverseDiff(true)),100; thinning=1, initial_params = paramvec)

    serialize("./fitted_models/samps_$(now()).dat",mhsamp)

    usdiagplots(alldata,mhsamp.value.data[50,1:end-1,1])
    usdiagplots(alldata,mhsamp.value.data[75,1:end-1,1])
    usdiagplots(alldata,mhsamp.value.data[100,1:end-1,1])

    hmcsamp = Turing.sample(alldata.model,HMC(1e-5,50; adtype = AutoReverseDiff(true)),20; init_params = paramvec)
#    nutsamp = Turing.sample(alldata.model,NUTS(300,.75; adtype=AutoReverseDiff(false)),100)
#    vfit = vi(alldata.model,ADVI(15,200),AutoReverseDiff(true))
end

end



function usdiagplots(alldata,paramvec,parnames)

    bar(1:length(paramvec),paramvec; title="Parameter Values",xlab="index") |> display
    #mapest=deserialize("fitted_models/USmodel_map_2024-08-29T17:24:15.073.dat")
    (densmin,densmax) = (alldata.model.args.densmin, alldata.model.args.densmax)
    (xmin,xmax) = (alldata.model.args.xmin, alldata.model.args.xmax)
    (ymin,ymax) = (alldata.model.args.ymin, alldata.model.args.ymax)

    kdindx = (1:36) .+ 6
    desindx = (1:36) .+ (6+36)
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),paramvec[kdindx] ./ 10)
    desirfun = Fun(ApproxFun.Chebyshev(xmin .. xmax) * ApproxFun.Chebyshev(ymin .. ymax ),paramvec[desindx] ./ 10)

    plotlyjs()# set the background

    p1=heatmap(kdfun;  title="dens Fun (dens,dens)",bins=(40,40),size=(600,600),xlim=(-7,7),ylim=(-7,7)) 

    p2=StatsPlots.histogram2d(StatsBase.sample(alldata.model.args.logreldens,1000),StatsBase.sample(alldata.model.args.logreldens,1000),normalize=true,
        size=(600,600),xlim=(-7,7),ylim=(-7,7),alpha=.25, title="Distribution of logreldens pairs")
    StatsPlots.plot(p1,p2; layout=(2,1)) |> display

    PlotlyJS.plot(PlotlyJS.histogram2dcontour(x=StatsBase.sample(alldata.model.args.logreldens,1000), y=StatsBase.sample(alldata.model.args.logreldens,1000)),
        PlotlyJS.Layout(alpha=0.1, title="Distribution of Randomly chosen Density Pairs",xlab="logreldensity",ylab="logreldensity")) |> display

    heatmap(desirfun, title="Desirability Fun (long,lat)",c=:rainbow) |> display

    preds = generated_quantities(alldata.model,paramvec,parnames)
    alldata.flows.preds = preds

    netactual = usnetmig(levelcode.(alldata.flows.fromcounty),levelcode.(alldata.flows.tocounty),alldata.flows.COUNT)
    netpred = usnetmig(levelcode.(alldata.flows.fromcounty),levelcode.(alldata.flows.tocounty),alldata.flows.preds)

    density(netactual-netpred; title="Residual Net error (actual - pred)",xlab="Count of people") |> display
    density((netactual-netpred) ./ alldata.geog.POPESTIMATE2016; xlim=(-.3,.3),legend=false, title="Residual Net error (actual - pred)",xlab="Fraction of source county") |> display
    density(netactual ./ alldata.geog.POPESTIMATE2016; xlim=(-.3,.3),legend=false,title="Net migration as fraction of population") |> display
    density(netactual ./ 1000.0; title="Actual Net Migration (thousands)",xlim=(-20,20)) |> display
    density(netactual ./ 1000.0; title="Actual Net Migration (thousands, full range") |> display

    scatter(netpred ./ 1000, netactual ./ 1000; ylab="actual net migration (thousands)",xlab="predicted (thousands)",xlim=(-20,20),ylim=(-20,20),
        title= "Net Migration prediction")
    Plots.abline!(1,0; legend=false) |> display

    scatter(netpred ./ 1000, netactual ./ 1000; ylab="actual net migration (thousands)",xlab="predicted (thousands)",
        title = "Net Migration prediction")
    Plots.abline!(1,0) |> display

    samps = StatsBase.sample(eachindex(alldata.flows.dist),5000; replace=false)

    flowsamp = alldata.flows[samps,:]

    density(log.(flowsamp.preds); label = "log(preds)", title="Density of log predictions")    
    density!(log.(flowsamp.COUNT .+ .01); label="actuals + .01") |> display
    StatsPlots.scatter(flowsamp.dist,log.((flowsamp.COUNT .+ .01) ./ flowsamp.preds); alpha=0.1, title="log(flow/pred) vs dist (in km)") |> display

    scatter(log.(flowsamp.preds/median(flowsamp.preds)),log.(flowsamp.COUNT/median(flowsamp.preds)); title="Comparing prediction to actual",
        xlab="log(pred/median(pred))", ylab="log(COUNT/median(pred))",xlim=(-5,5),ylim=(-5,5), alpha=0.1) |> display


    StatsPlots.scatter(flowsamp.dist,[log.((r.COUNT .+ 1) ./alldata.geog.POPESTIMATE2016[levelcode(r.fromcounty)]) for r in eachrow(flowsamp)];
        title = "log((flow+1)/from_pop) vs distance", alpha=0.1) |> display

    StatsPlots.scatter(flowsamp.dist,[log.((r.preds) ./alldata.geog.POPESTIMATE2016[levelcode(r.fromcounty)]) for r in eachrow(flowsamp)];
        title = "log((pred+1)/from_pop) vs distance", alpha=0.1) |> display

    kfrom = paramvec[6] / 10.0
    density(kfrom .* log.(alldata.geog.POPESTIMATE2016[levelcode.(flowsamp.tocounty)] ./ median(alldata.geog.POPESTIMATE2016)); title="DIstribution of log population US * kfrom") |> display

    tots = totalflow(alldata.flows.fromcounty,alldata.flows.tocounty,alldata.flows.COUNT,alldata.geog.countyid)
    totpreds = totalflow(alldata.flows.fromcounty,alldata.flows.tocounty,alldata.flows.preds,alldata.geog.countyid)

    scatter(log.(totpreds ./ alldata.geog.POPESTIMATE2016),log.(tots ./ alldata.geog.POPESTIMATE2016),
        xlab="log(Pred total flux / Pop)", ylab = "log(Actual Total Flux / Pop)", title = "Total Flux comparison")
    Plots.abline!(1,0; label = "y = x", legend = false) |> display

end

grabparams(chain,n) = chain.value.data[n,1:end-1,1]


function fitandwritefile(alldata,flowout,geogout,densout,paramout)

    algo = LBFGS()

    parnames = [["a","c","d0","dscale","neterr","ktopop"];
        ["kd[$i]" for i in 1:36];
        ["desirecoefs[$i]" for i in 1:36]]


    lb = [[-60.0,0.0,0.0,1.0,0.01,-10.0];
            fill(-50.50,36);
            fill(-50.50,36)]
    ub = [[60.0, 20.0,10.0,15.0,100.0,10.0];
            fill(50.50,36);
            fill(50.50,36)]
    ini = rand(Normal(0.0,0.10),length(ub))
    ini[1:7] .= [-7.6,1.81,1.5,5.0,1.0,3.5,0.0]
    println("Optimization starts")
    mapest = maximum_a_posteriori(alldata.model, algo; adtype = AutoReverseDiff(false),
                                  initial_params = ini, lb = lb, ub = ub, maxiters = 500, maxtime = 400,
                                  reltol=1e-3, progress=true)
    println("Optimization finished")
    paramvec = mapest.values.array

    #usdiagplots(alldata,paramvec,parnames)
    println("Sampling starts")
    mhsamp = Turing.sample(alldata.model, MH(.1^2*I(length(mapest.values))), 10;
                           thinning=10, initial_params = paramvec)
    println("Sampling finished")
#    mhsamp = Turing.sample(alldata.model,HMCDA(200,.7,1.0; adtype=AutoReverseDiff(true)),100; thinning=1, initial_params = paramvec)
    paramvec = grabparams(mhsamp,10)

    preds = generated_quantities(alldata.model,paramvec,parnames)
    alldata.flows.preds = preds

    CSV.write(flowout,alldata.flows)

    (densmin,densmax) = (alldata.model.args.densmin, alldata.model.args.densmax)
    (xmin,xmax) = (alldata.model.args.xmin, alldata.model.args.xmax)
    (ymin,ymax) = (alldata.model.args.ymin, alldata.model.args.ymax)    
    kdindx = (1:36) .+ 6
    desindx = (1:36) .+ (6+36)
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),paramvec[kdindx] ./ 10)
    desirfun = Fun(ApproxFun.Chebyshev(xmin .. xmax) * ApproxFun.Chebyshev(ymin .. ymax ),paramvec[desindx] ./ 10)

    alldata.geog.desirability = [desirfun(x,y) for (x,y) in zip(alldata.geog.x,alldata.geog.y)]

    CSV.write(geogout,alldata.geog)
    densvals = range(minimum(alldata.geog.logreldens),maximum(alldata.geog.logreldens),100)
    densfundf = DataFrame((fromdens=fd, todens=td, funval=kdfun(fd,td)) for fd in densvals, td in densvals)
    CSV.write(densout,densfundf)
    CSV.write(paramout,DataFrame(paramval=paramvec,parname=parnames))
end


function main()
    @sync begin
        Threads.@spawn begin
            usd = loadallUSdata(0; sample = sample) # add no zeros, 

            Random.seed!(Int64(datetime2unix(DateTime("2024-10-01T09:07:14")))) # seed based on current time when I wrote the function
            usd.geog.x = usd.geog.INTPTLONG
            usd.geog.y = usd.geog.INTPTLAT
            fitandwritefile(usd,"manuscript_input/usflows.csv","manuscript_input/usgeog.csv","manuscript_input/usdensfun.csv","manuscript_input/usparams.csv")
        end
        germ = loadallGermData(0; sample = sample)
        germ.geog.x = germ.geog.xcoord
        germ.geog.y = germ.geog.ycoord

        Threads.@threads for age in levels(germ.flows.agegroup)
            agedat = @subset(germ.flows,germ.flows.agegroup .== age)
            modl = usmodel(agedat.flows,sum(agedat.flows),levelcode.(agedat.fromdist), levelcode.(agedat.todist),
                            median(germ.geog.pop),agedat.dist,
                            germ.geog.xcoord,minimum(germ.geog.xcoord),maximum(germ.geog.xcoord),
                            germ.geog.ycoord,minimum(germ.geog.ycoord),maximum(germ.geog.ycoord),
                            germ.geog.logreldens,minimum(germ.geog.logreldens),maximum(germ.geog.logreldens),
                            germ.geog.pop,nrow(germ.geog),100.0,36,36) ## nothing == netactual, we're not using it anymore
            fitandwritefile((flows=agedat,geog=germ.geog,model=modl),
                            "manuscript_input/germflows_$(age).csv",
                            "manuscript_input/germgeog_$(age).csv",
                            "manuscript_input/germdensfun_$(age).csv",
                            "manuscript_input/germparams_$(age).csv")
        end
    end
    println("Computation finished!")
end


# getUSflows()
# getUSgeog()
# getUScountypop()
sample = false
main()
