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
net = calc_net2(df, :flows, [:agegroup])
leftjoin!(net, year(di, 2017)[!, [:distcode, :xcoord, :ycoord]],
          on = :fromdist => :distcode)

function nmrage(f, net, shp, st, a, x, y)
    df = age(net, a)
    ax = Axis(f[x, y], aspect=DataAspect(),
          title = "$a")
    viz!(ax, rmreg(shp, :AGS).geometry, color = df.nmr, colorrange = extrema(net.nmr));
    hideall!(ax)
    overlay_states(ax, st)
end

f = Figure(size = (700, 700), fontsize = 12);
nmrage(f, net, shp, st, "below18", 1, 1);
nmrage(f, net, shp, st, "18-25", 1, 2);
nmrage(f, net, shp, st, "25-30", 1, 3);
nmrage(f, net, shp, st, "30-50", 2, 1);
nmrage(f, net, shp, st, "50-65", 2, 2);
nmrage(f, net, shp, st, "above65", 2, 3);
Colorbar(f[3, 1:3], limits = extrema(net.nmr), vertical = false,
         height = 3, label = "Net Migration Rate", width = Relative(.9))

titlelayout = GridLayout(f[0, 1], tellwidth = false)
Label(titlelayout[1, 1], "Net Migration Rate", fontsize = 15, font = "TeX Gyre Heros Bold Makie")
Label(titlelayout[2, 1], "Average 2000-2017", fontsize = 10)
rowgap!(titlelayout, 0)
rowgap!(f.layout, 0)
colgap!(f.layout, -50)
resize_to_layout!(f)
save(joinpath(outp, "nmr_age.pdf"), f)
