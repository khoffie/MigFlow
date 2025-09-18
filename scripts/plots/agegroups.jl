using DataFrames, CSV, CairoMakie, Statistics, StatsBase, Revise
using LaTeXStrings, AlgebraOfGraphics, GeoIO, GeoStats, GLM
using KernelDensity
outp = "/home/konstantin/paper/sections/texdocs/images/"

includet("../../src/diagplots.jl")
includet("../../src/analyzegeo.jl")
includet("../../src/model_helpers.jl")
includet("../../src/samplecircle.jl")

df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
df = rmreg(df, :fromdist)
df = rmreg(df, :todist)

di = CSV.read("../data/districts.csv", DataFrame)
di.area = di.pop ./ di.density;

pop = combine(groupby(unique(df, [:fromdist, :year, :agegroup]),
                      [:year, :agegroup]), :frompop => sum => :agepop)
## combine(groupby(pop, :year), :agepop => sum)
shp = GeoIO.load("../../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp")

years = unique(df.year)
ages = unique(df.agegroup)


######################################################################
################# total and relative flows ###########################
dft = combine(groupby(df, [:year, :agegroup]), :flows => sum => :total)
leftjoin!(dft, pop, on = [:year, :agegroup])
dft.prob = dft.total ./ dft.agepop .* 100
dft.total = dft.total ./ 100e3

ticks = [2000, 2005, 2010, 2015, 2017]
vs = ["2000", "05", "10", "15", "17"]

f = Figure(size = (300, 150), fontsize = 10);
ax1 = Axis(f[1, 1],
           ylabel = "Flows (100k)",
           xlabel = "Year",
           xticks = (ticks, vs),
           xgridvisible = false, ygridvisible = false);
ax2 = Axis(f[1, 2],
           ylabel = "Movement Frequency (%)",
           xlabel = "Year",
           xticks = (ticks, vs),
           xgridvisible = false, ygridvisible = false);
for a in ages
    foo = age(dft, [a])
    lines!(ax1, foo.year, foo.total, label = a)
    lines!(ax2, foo.year, foo.prob, label = a)

end

Legend(f[1, 3], ax1, patchsize = (2, 10), framevisible = false, labelsize = 10)
resize_to_layout!(f)
save(joinpath(outp, "flows_age.pdf"), f)
