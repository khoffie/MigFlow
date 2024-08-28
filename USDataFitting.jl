using CSV,DataFrames,StatsPlots,Distributions,Turing,StatsBase,StatsFuns,FixedWidthTables,DataFramesMeta
using Printf, CategoricalArrays
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

const non48fips = (2,11,15,72) #fips for AK, HI, DC, PR

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



function loadallUSdata()

    flows = loadUS48flows()
    geog = loadUS48geog()
    pop = loadUScountypop()
    pop = DataFramesMeta.@select!(pop,:STATEFP,:COUNTYFP,:POPESTIMATE2016)
    geog = leftjoin(geog,pop,on = [:STATEFP,:COUNTYFP])

    meddens = median(geog.POPESTIMATE2016 ./ geog.ALAND)
    geog.logreldens = log.(geog.POPESTIMATE2016 ./ geog.ALAND / meddens)
    geog

    flowset = Set((f.fromcounty,f.tocounty) for f in eachrow(flows))
    zerosamp = DataFrame((COUNT = 0,fromcounty = f,tocounty = t) for (f,t) in 
        zip(sample(geog.countyid,10000),sample(geog.countyid,10000)) if !((f,t)  in flowset))

    (geog=geog, county = geog, flows = flows, zerosamp=zerosamp)
end

#=

densdict = Dict(alldata.county.countyid .=> alldata.county.logreldens)
histogram2d([densdict[f] for f in alldata.flows.fromcounty],[densdict[f] for f in alldata.flows.tocounty])

=#