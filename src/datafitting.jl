# using CSV, DataFrames, FixedWidthTables, DataFramesMeta, CategoricalArrays
# using LibGit2
# using StatsBase, StatsFuns, StatsPlots, Distributions, Random, LaTeXStrings, Plots
# using Turing, ReverseDiff, ApproxFun
# using Printf, Revise, Dates, Enzyme, Serialization, SliceSampling
# using LogDensityProblems, LogDensityProblemsAD, Distances, LinearAlgebra
# Enzyme.API.runtimeActivity!(true) ## to deal with an Enzyme bug, per https://discourse.julialang.org/t/enzyme-ready-for-everyday-use-2024/118819/7
# import PlotlyJS

# includet("models.jl")
# includet("samplerows.jl")
# includet("fitandwrite.jl")
# includet("temperedmodel.jl")
# includet("utils.jl")
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

    if settings[:rm_dir] && isdir(outpath)
        rm(outpath, recursive = true)
    end
    mkpath(outpath)
    savesettings(settings, outpath)
    file = "./writeup/juliaout_path.txt" ## report.Rmd reads the path from here
    write(file, outpath)

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
        
        germ = loadallGermData(0; sample = settings[:sample_rows],
                               positive_only = settings[:positive_only])
        germ.geog.x = germ.geog.xcoord
        germ.geog.y = germ.geog.ycoord

        ##        Threads.@threads for age in levels(germ.flows.agegroup)
        ages = settings[:agegroups] == nothing ? levels(germ.flows.agegroup) : settings[:agegroups]
        Threads.@threads for age in ages
            for year in unique(germ.flows.year)
                agedat = @subset(germ.flows, :agegroup .== age, :year .== year)
                geodat = @subset(germ.geog, :year .== year)
                mdl = germmodel(agedat, geodat,
                                settings[:model_type],
                                settings[:positive_only])
                outpaths = createpaths(outpath, "germ", year, age)
                germd = (flows = agedat, geog = geodat, model = mdl)
                settings[:fit_germ] ? fitandwritefile(germd, settings, outpaths) : println("German data not fitted")
            end
        end
    end
    println("Computation finished!")
end
# getUSflows()
# getUSgeog()
# getUScountypop()

function makeoutpath(outpath)
    datepath = joinpath("results", Dates.format(now(), "yyyy-mm-dd_HH-MM-SS"))
    out = isnothing(outpath) ? datepath : joinpath("results", outpath)
    println("Output saved into " * out)
    return out
end
