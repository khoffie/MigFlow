using Pkg
Pkg.activate(".")
using CSV, DataFrames, FixedWidthTables, DataFramesMeta, CategoricalArrays
using LibGit2
using StatsBase, StatsFuns, StatsPlots, Distributions, Random, LaTeXStrings, Plots
using Turing, ReverseDiff, ApproxFun
using Printf, Revise, Dates, Enzyme, Serialization, SliceSampling
using LogDensityProblems, LogDensityProblemsAD, Distances, LinearAlgebra
Enzyme.API.runtimeActivity!(true) ## to deal with an Enzyme bug, per https://discourse.julialang.org/t/enzyme-ready-for-everyday-use-2024/118819/7
import PlotlyJS

includet("models.jl")
includet("samplerows.jl")
includet("postprocess.jl")
includet("fitandwrite.jl")
includet("temperedmodel.jl")
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

function loadGermFlows(distcodes, positive_only)
    fl = CSV.read("data/FlowDataGermans.csv", DataFrame)
    fl = filter(row -> row.fromdist in distcodes, fl)
    fl = filter(row -> row.todist in distcodes, fl)
    fl.fromdist = categorical(fl.fromdist, levels = distcodes)
    fl.todist = categorical(fl.todist, levels = distcodes)
    if positive_only
        fl = filter(row -> row.flows .> 0, fl)
    end
    return fl
end


function loadallGermData(nzeros = 0; sample, positive_only)
    geog = loadGermGeog()
    if sample == true
        geog = sample_germ(geog)
    end
    flows = loadGermFlows(levels(geog.distcode), positive_only)
    return (geog=geog, flows=flows)
end

function loadallUSdata(nzeros = 0; sample, positive_only)

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
    # if positive_only
    #     flows = filter(row -> row.flows .> o, flows)        
    # end
    flows.dist = [haversine((geog[levelcode(r.fromcounty),:INTPTLONG],geog[levelcode(r.fromcounty),:INTPTLAT]), 
                        (geog[levelcode(r.tocounty),:INTPTLONG],geog[levelcode(r.tocounty),:INTPTLAT])) / 1000.0 for r in eachrow(flows)] ## haversine is in m, calculate km
    #netactual::Vector{Float64} = usnetmig(levelcode.(flows.fromcounty),levelcode.(flows.tocounty),flows.COUNT)

    mod = usmodel(flows.COUNT,sum(flows.COUNT),levelcode.(flows.fromcounty),
                levelcode.(flows.tocounty),median(geog.POPESTIMATE2016),flows.dist,
                geog.INTPTLONG,minimum(geog.INTPTLONG),maximum(geog.INTPTLONG),
                geog.INTPTLAT,minimum(geog.INTPTLAT),maximum(geog.INTPTLAT),
                geog.logreldens,minimum(geog.logreldens),maximum(geog.logreldens),
                geog.POPESTIMATE2016,nrow(geog),100.0,36,36, positive_only)
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

function createpaths(path, type, year, age)
    datasets = ["flows", "geog", "densfun", "params", "chain"]
    paths = Dict(d => "$(path)/$(type)$(d)_$(year)_$(age).csv" for d in datasets)
    return paths
end

function mainfit(settings, outpath)
    function savesettings(settings, path)
        settings = DataFrame(setting = [string(k) for k in keys(settings)],
                             value = [string(k) for k in values(settings)])
        CSV.write(joinpath(path, "settings.csv"), settings)
    end

    mkpath(outpath)
    @sync begin
        Threads.@spawn begin
            usd = loadallUSdata(0; sample = settings[:sample_rows],
                                positive_only = settings[:positive_only]) # add no zeros, 

            Random.seed!(Int64(datetime2unix(DateTime("2024-10-01T09:07:14")))) # seed based on current time when I wrote the function
            usd.geog.x = usd.geog.INTPTLONG
            usd.geog.y = usd.geog.INTPTLAT
            outpaths = createpaths(outpath, "us", 999,"all")
            settings[:fit_us] ? fitandwritefile(usd, settings, outpaths) : println("US data not fitted")
        end
        
        germ = loadallGermData(0; sample = settings[:sample_rows], positive_only = settings[:positive_only])
        germ.geog.x = germ.geog.xcoord
        germ.geog.y = germ.geog.ycoord

        ##        Threads.@threads for age in levels(germ.flows.agegroup)
        ages = settings[:agegroups] == nothing ? levels(germ.flows.agegroup) : settings[:agegroups]
        Threads.@threads for age in ages
            for year in unique(germ.flows.year)
                agedat = @subset(germ.flows, :agegroup .== age, :year .== year)
                geodat = @subset(germ.geog, :year .== year)
                mdl = germmodel(agedat, geodat, settings[:positive_only])
                outpaths = createpaths(outpath, "germ", year, age)
                germd = (flows = agedat, geog = geodat, model = mdl)
                settings[:fit_germ] ? fitandwritefile(germd, settings, outpaths) : println("German data not fitted")
            end
        end
    end
    println("Computation finished!")
    savesettings(settings, outpath)

    file = "./writeup/juliaout_path.txt" ## report.Rmd reads the path from here
    open(file, "w") do f
        write(f, outpath)
   end
end


# getUSflows()
# getUSgeog()
# getUScountypop()

function makeoutpath(outpath)
    datepath = joinpath("manuscript_input", Dates.format(now(), "yyyy-mm-dd_HH-MM-SS"))
    out = isnothing(outpath) ? datepath : joinpath("manuscript_input", outpath)
    println("Output saved into " * out)
    return out
end
