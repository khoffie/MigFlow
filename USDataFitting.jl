using Pkg
Pkg.activate(".")
using CSV, DataFrames, FixedWidthTables, DataFramesMeta, CategoricalArrays, RCall
using StatsBase, StatsFuns, StatsPlots, Distributions, Random, StatProfilerHTML
using Turing, OptimizationOptimJL, ApproxFun, OptimizationBBO, OptimizationNLopt, NLopt, ReverseDiff
using Printf, Revise, Dates, Enzyme, Serialization
using LogDensityProblems, LogDensityProblemsAD, Distances, LinearAlgebra
Enzyme.API.runtimeActivity!(true) ## to deal with an Enzyme bug, per https://discourse.julialang.org/t/enzyme-ready-for-everyday-use-2024/118819/7
import PlotlyJS

includet("models.jl")
includet("samplerows.jl")
includet("post_process.jl")
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


grabparams(chain,n) = chain.value.data[n,1:end-1,1]

function fitandwritefile(alldata, flowout, geogout, densout, paramout, chainout)
    function gen_inits()
        parnames = [["a", "c", "d0", "dscale", "ktopop"];
                    ["kd[$i]" for i in 1:36];
                    ["desirecoefs[$i]" for i in 1:36]]
        lb = [[-60.0, 0.0, 0.0, 1.0, -10.0];
              fill(-50.50, 36);
              fill(-50.50, 36)]
        ub = [[60.0, 20.0, 10.0, 15.0, 10.0];
              fill(50.50, 36);
              fill(50.50, 36)]
        ini = rand(Normal(0.0, 0.10), length(ub))
        ini[1:6] .= [-7.6, 1.81, 1.5, 5.0, 3.5, 0.0]
        df = DataFrame(:pars => parnames, :lb => lb, :ub => ub, :inits => ini)
        return(df)
    end

    function runoptim(vals, optim)
        if optim == true
            algo = LBFGS()
            println("Optimization starts")
            mapest = maximum_a_posteriori(alldata.model, algo; adtype = AutoReverseDiff(false),
                                          initial_params = vals.inits, lb = vals.lb, ub = vals.ub,
                                          maxiters = 50, maxtime = 400,
                                          reltol=1e-3, progress = true)
            println("Optimization finished")
            vals.optis = mapest.values.array
        else
            println("No optimization, using random inits for sampling")
            vals.optis = vals.inits
        end
        println(vals[[1:10; 43:47], :])
        return vals
    end

    vals = gen_inits()
    vals = runoptim(vals, false)
    #usdiagplots(alldata,paramvec,parnames)
    println("Sampling starts")
    mhsamp = Turing.sample(alldata.model, MH(.1^2*I(length(mapest.values))), 10;
                           thinning = 5, initial_params = vals.optis)
    # plt = Plots.plot(mhsamp[:lp])
    # Plots.savefig(plt, plotout)

    Serialization.serialize(chainout, mhsamp)
    println("Sampling finished")
#    mhsamp = Turing.sample(alldata.model,HMCDA(200,.7,1.0; adtype=AutoReverseDiff(true)),100; thinning=1, initial_params = paramvec)
    vals.optsam = grabparams(mhsamp,10)
    println(vals[[1:10; 43:47], :])
    preds = generated_quantities(alldata.model, vals.optsam, vals.pars)
    alldata.flows.preds = preds

    CSV.write(flowout,alldata.flows)

    (densmin,densmax) = (alldata.model.args.densmin, alldata.model.args.densmax)
    (xmin,xmax) = (alldata.model.args.xmin, alldata.model.args.xmax)
    (ymin,ymax) = (alldata.model.args.ymin, alldata.model.args.ymax)    
    kdindx = (1:36) .+ 5 #number of non kd or cheby inits/ priors
    desindx = (1:36) .+ (5 + 36)
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
                vals.optsam[kdindx] ./ 10)
    desirfun = Fun(ApproxFun.Chebyshev(xmin .. xmax) * ApproxFun.Chebyshev(ymin .. ymax ),
                   vals.optsam[desindx] ./ 10)

    alldata.geog.desirability = [desirfun(x,y) for (x,y) in zip(alldata.geog.x,alldata.geog.y)]

    CSV.write(geogout,alldata.geog)
    densvals = range(minimum(alldata.geog.logreldens),maximum(alldata.geog.logreldens),100)
    densfundf = DataFrame((fromdens=fd, todens=td, funval=kdfun(fd,td)) for fd in densvals, td in densvals)
    CSV.write(densout,densfundf)
    CSV.write(paramout,DataFrame(paramval = vals.optsam, parname = vals.pars))
end


function main()
    @sync begin
        Threads.@spawn begin
            usd = loadallUSdata(0; sample = sample) # add no zeros, 

            Random.seed!(Int64(datetime2unix(DateTime("2024-10-01T09:07:14")))) # seed based on current time when I wrote the function
            usd.geog.x = usd.geog.INTPTLONG
            usd.geog.y = usd.geog.INTPTLAT
            fitandwritefile(usd,
                            "manuscript_input/usflows.csv",
                            "manuscript_input/usgeog.csv",
                            "manuscript_input/usdensfun.csv",
                            "manuscript_input/usparams.csv",
                            "manuscript_input/uschain")
        end
        germ = loadallGermData(0; sample = sample)
        germ.geog.x = germ.geog.xcoord
        germ.geog.y = germ.geog.ycoord

        Threads.@threads for age in levels(germ.flows.agegroup)
            println(age)
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
                            "manuscript_input/germparams_$(age).csv",
                            "manuscript_input/germchain_$(age)")
        end
    end
    println("Computation finished!")
end


# getUSflows()
# getUSgeog()
# getUScountypop()
sample = false
main()
post_process()
## R"helpeR::render_doc('~/Documents/GermanMigration/writeup', 'report.Rmd')"
## R"helpeR::render_doc('~/Documents/GermanMigration/writeup', 'math.tex')"
