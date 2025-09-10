using CSV, DataFrames, Revise, GeoStats, CairoMakie, StatsBase
includet("../src/analyzegeo.jl")
includet("plots/utils.jl")
function agelines!(df, ax, ages, xcol, ycol)
    for a in ages
        foo = age(df, a)
        lines!(ax, foo[!, xcol], foo[!, ycol], label = a)
    end
end


df = CSV.read("../data/FlowDataGermans.csv", DataFrame);
districts = CSV.read("../data/districts.csv", DataFrame);
correct = CSV.read("../data/clean/correct.csv", DataFrame)
correct = year(correct, years)
ages = unique(df.agegroup)
allyears = unique(df.year)

## In the early 2000s are unreasonably high flows. What is going on?
## To get a quick idea we will count the flows that are 5x higher than
## the average between 2006-2008

df2 = outflux(df, sum, [:todist, :year])
dub = df2[df2.year .<= 2005, :]
ave = outflux(year(df2, 2006.0, 2008.0), mean, [:todist], "average")
leftjoin!(dub, ave, on = [:fromdist, :todist])
dub.rel = dub.flows ./ dub.average
dub.toolarge = dub.rel .>= 5.0
df3 = sort(combine(groupby(dub, [:fromdist, :year]), :toolarge => sum => :N), :N)
sort(combine(groupby(df3, :fromdist), :N => sum => :N), :N)


f = Figure(size = (600, 400), fontsize = 10);
ax1 = Axis(f[1, 1],
           title = "Outflux district Göttingen",
           xlabel = "Year",
           ylabel = "Outflux (in thousands)")
lines!(ax1, allyears, origin(outflux(df, sum, [:year]), 3159).flows ./ 1000)
ax2 = Axis(f[1, 2],
           title = "Flows from Göttingen to Unna",
           xlabel = "Year",
           ylabel = "Flows (in thousands)")
lines!(ax2, allyears, destination(origin(df2, :3159), :5978).flows ./ 1000)
ax3 = Axis(f[1, 3],
           title = "Flows from Göttingen to Spree-Neiße",
           xlabel = "Year",
           ylabel = "Flows (in thousands)")
lines!(ax3, allyears, destination(origin(df2, :3159), :12071).flows ./ 1000)


influx(destination(df2, 3159), sum, [:year])
influx(destination(df2, 5978), sum, [:year])

code(districts, 12071)
code(districts, 5978)

filter(row -> row.ags_new in x, correct)
getname(df, x) = unique(df[df.ags_new .== x, :name_old])[1]
