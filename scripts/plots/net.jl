using DataFrames, CSV, CairoMakie, Statistics, StatsBase, Revise
using LaTeXStrings, AlgebraOfGraphics, GeoIO, GeoStats, GLM
using KernelDensity
outp = "/home/konstantin/paper/sections/texdocs/images/"

includet("../src/diagplots.jl")
includet("../src/analyzegeo.jl")
includet("../src/model_helpers.jl")
includet("../src/samplecircle.jl")

df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
df = rmreg(df, :fromdist)
df = rmreg(df, :todist)

di = CSV.read("../data/districts.csv", DataFrame)
di.area = di.pop ./ di.density;

pop = combine(groupby(unique(df, [:fromdist, :year, :agegroup]),
                      [:year, :agegroup]), :frompop => sum => :agepop)
## combine(groupby(pop, :year), :agepop => sum)
shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp")

years = unique(df.year)
ages = unique(df.agegroup)

######################################################################
net = calc_net(df, :flows)
sort(net, :nmr)

function calc_net2(df, col, group)
    netf = combine(DataFrames.groupby(df, [:fromdist, group...]), col => sum)
    rename!(netf, string(col) * "_sum" => :outflux)
    nett = combine(DataFrames.groupby(df, [:todist, group...]), col => sum)
    rename!(nett, string(col) * "_sum" => :influx)
    net = innerjoin(netf, nett, on = [:fromdist => :todist, group...])
    net.net = net.influx .- net.outflux
    net.total = net.influx .+ net.outflux
    net.nmr = net.net ./ net.total
    return net
end
net = calc_net2(df, :flows, [:agegroup, :year])
sort(net, :nmr)
origin(net, [7211])

sum(origin(age(year(df, 2006), "18-25"), [7211]).flows)
sort(destination(age(year(df, 2006), "18-25"), [7211]), :flows)

sort(calc_net(df, :flows), :nmr)
