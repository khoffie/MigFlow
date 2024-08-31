using CSV,DataFrames,StatsPlots,Distributions,Turing,StatsBase,StatsFuns,FixedWidthTables,DataFramesMeta

using Enzyme
Enzyme.API.runtimeActivity!(true) ## to deal with an Enzyme bug, per https://discourse.julialang.org/t/enzyme-ready-for-everyday-use-2024/118819/7

using Printf, CategoricalArrays
import PlotlyJS
using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random,
    ReverseDiff, Revise, RCall, Dates
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta,
    StatProfilerHTML, StatsFuns, OptimizationBBO, Printf,OptimizationNLopt,NLopt, LinearAlgebra

using LogDensityProblems,LogDensityProblemsAD, Distances

includet("models.jl")
include("models.jl")

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



function loadallUSdata(nzeros = 100000)

    flows = loadUS48flows()
    geog = loadUS48geog()
    pop = loadUScountypop()
    pop = DataFramesMeta.@select!(pop,:STATEFP,:COUNTYFP,:POPESTIMATE2016)
    geog = leftjoin(geog,pop,on = [:STATEFP,:COUNTYFP])
    meddens = median(geog.POPESTIMATE2016 ./ geog.ALAND)
    geog.logreldens = log.(geog.POPESTIMATE2016 ./ geog.ALAND / meddens)
    geog.countylevel = levelcode.(geog.countyid)
    sort!(geog,[:countylevel])
    levels!(flows.fromcounty,levels(geog.countyid))
    levels!(flows.tocounty,levels(geog.countyid))

#    flowset = Set((f.fromcounty,f.tocounty) for f in eachrow(flows))
#    zerosamp = DataFrame((COUNT = 0,fromcounty = f,tocounty = t) for (f,t) in 
#        zip(StatsBase.sample(geog.countyid,nzeros),StatsBase.sample(geog.countyid,nzeros)) if !((f,t)  in flowset))
#    flows = [flows; zerosamp]

    flows.dist = [haversine((geog[levelcode(r.fromcounty),:INTPTLONG],geog[levelcode(r.fromcounty),:INTPTLAT]), 
                        (geog[levelcode(r.tocounty),:INTPTLONG],geog[levelcode(r.tocounty),:INTPTLAT])) / 1000.0 for r in eachrow(flows)] ## haversine is in m, calculate km
    mod = usmodel(flows.COUNT,sum(flows.COUNT),levelcode.(flows.fromcounty),
                levelcode.(flows.tocounty),median(geog.POPESTIMATE2016),flows.dist,
                geog.INTPTLONG,minimum(geog.INTPTLONG),maximum(geog.INTPTLONG),
                geog.INTPTLAT,minimum(geog.INTPTLAT),maximum(geog.INTPTLAT),
                geog.logreldens,minimum(geog.logreldens),maximum(geog.logreldens),
                geog.POPESTIMATE2016,nrow(geog),100.0,[],36,36)
    (geog=geog, flows = flows, model=mod) #zerosamp=zerosamp,model=mod)
end

if false
    densdict = Dict(alldata.county.countyid .=> alldata.county.logreldens)
    histogram2d([densdict[f] for f in alldata.flows.fromcounty],[densdict[f] for f in alldata.flows.tocounty])


    alldata = loadallUSdata()

    #algo = BBO_de_rand_1_bin_radiuslimited()
    #algo = LBFGS()
    algo = NLopt.LN_BOBYQA()

    parnames = [["a","c","d0","dscale","neterr","ktopop"];
        ["kd[$i]" for i in 1:36];
        ["desirecoefs[$i]" for i in 1:36]]


    lb = [[-30.0,0.0,0.0,1.0,0.01,-10.0];
            fill(-20.50,36);
            fill(-20.50,36)]
    ub = [[30.0, 5.0,10.0,15.0,4.0,10.0];
            fill(20.50,36);
            fill(20.50,36)]
    ini = rand(Normal(0.0,0.10),length(ub))
    ini[1:7] .= [-9.6,1.81,1.5,5.0,1.0,3.5,0.0] 
    mapest = maximum_a_posteriori(alldata.model, algo; adtype = AutoReverseDiff(),initial_params=ini,
        lb=lb,ub=ub, maxiters = 800, maxtime = 600, reltol=1e-3,progress=true)
    serialize("fitted_models/USmodel_map_$(now()).dat",mapest)

    
    bar(1:length(mapest.values.array),mapest.values.array; title="Parameter Values",xlab="index") |> display
    #mapest=deserialize("fitted_models/USmodel_map_2024-08-29T17:24:15.073.dat")
    (densmin,densmax) = (alldata.model.args.densmin, alldata.model.args.densmax)
    (xmin,xmax) = (alldata.model.args.xmin, alldata.model.args.xmax)
    (ymin,ymax) = (alldata.model.args.ymin, alldata.model.args.ymax)

    kdindx = (1:36) .+ 6
    desindx = (1:36) .+ (6+36)
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),mapest.values.array[kdindx] ./ 10)
    desirfun = Fun(ApproxFun.Chebyshev(xmin .. xmax) * ApproxFun.Chebyshev(ymin .. ymax ),mapest.values.array[desindx] ./ 10)

    plotlyjs()# set the background

    p1=heatmap(kdfun;  title="dens Fun (dens,dens)",bins=(40,40),size=(600,600),xlim=(-7,7),ylim=(-7,7)) 

    p2=StatsPlots.histogram2d(StatsBase.sample(alldata.model.args.logreldens,1000),StatsBase.sample(alldata.model.args.logreldens,1000),normalize=true,
        size=(600,600),xlim=(-7,7),ylim=(-7,7),alpha=.25, title="Distribution of logreldens pairs")
    StatsPlots.plot(p1,p2; layout=(2,1)) |> display

    PlotlyJS.plot(PlotlyJS.histogram2dcontour(x=StatsBase.sample(alldata.model.args.logreldens,1000), y=StatsBase.sample(alldata.model.args.logreldens,1000)),
        PlotlyJS.Layout(alpha=0.1, title="Distribution of Randomly chosen Density Pairs",xlab="logreldensity",ylab="logreldensity")) |> display

    heatmap(desirfun, title="Desirability Fun (long,lat)",c=:rainbow) |> display

    preds = generated_quantities(alldata.model,mapest.values.array,parnames)
    alldata.flows.preds = preds
    samps = StatsBase.sample(eachindex(alldata.flows.dist),5000; replace=false)

    flowsamp = alldata.flows[samps,:]

    density(log.(flowsamp.preds); label = "log(preds)", title="Density of log predictions")    
    density!(log.(flowsamp.COUNT .+ .01); label="actuals + .01") |> display
    StatsPlots.scatter(flowsamp.dist,log.(flowsamp.COUNT ./ flowsamp.preds); alpha=0.1, title="log(flow/pred) vs dist (in km)") |> display


    StatsPlots.scatter(flowsamp.dist,[log.((r.COUNT .+ 1) ./alldata.geog.POPESTIMATE2016[levelcode(r.fromcounty)]) for r in eachrow(flowsamp)];
        title = "log((flow+1)/from_pop) vs distance", alpha=0.1) |> display

    StatsPlots.scatter(flowsamp.dist,[log.((r.preds) ./alldata.geog.POPESTIMATE2016[levelcode(r.fromcounty)]) for r in eachrow(flowsamp)];
        title = "log((pred+1)/from_pop) vs distance", alpha=0.1) |> display

    kfrom = mapest.values.array[6] / 10.0
    density(kfrom .* log.(alldata.geog.POPESTIMATE2016[levelcode.(flowsamp.tocounty)] ./ median(alldata.geog.POPESTIMATE2016)); title="DIstribution of log population US * kfrom") |> display


#    vfit = vi(alldata.model,ADVI(15,200),AutoReverseDiff(true))
end
