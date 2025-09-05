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
####################### Flow Magnitude ###############################
dfm = sort(outflux(df, sum, [:todist]), :outflux)
dfm.yearly_flows = dfm.outflux ./ length(years)
dfdist = year(age(df, "18-25"), 2017)[!, [:fromdist, :todist, :dist]]
leftjoin!(dfm, dfdist, on = [:fromdist, :todist])
flows = dfm.yearly_flows

qs = quantile(flows, range(0, 1, 100))
flows_sum = [sum(flows[flows .> qs[i-1] .&& flows .< qs[i]]) for i in 2:length(qs)]
flows_sum = flows_sum ./ sum(flows) .* 100
df1 = DataFrame(qs = 1:(length(qs) - 1), flows = flows_sum)
ticks = [.25, .5, .75, .9, .95, 1] * 100
vs = Int.(round.(quantile(flows, ticks ./ 100), digits = 0))
vs = string.(Int.(ticks)) .* "\n" .* "(" .* string.(vs) .* ")"

######################################################################
############### ecdf most important regions ##########################
dfsout = combine(groupby(df, [:fromdist, :todist]),
              [:flows => sum => :flows])
leftjoin!(dfsout, outflux(df, sum), on = :fromdist)
dfsout.cond = dfsout.flows ./ dfsout.outflux .* 100
dfsout = sort(dfsout, [:fromdist, :cond], rev = true)
dfsout.id = repeat(1:399, 400)

dfsin = combine(groupby(df, [:todist, :fromdist]),
              [:flows => sum => :flows])
leftjoin!(dfsin, influx(df, sum), on = :todist)
dfsin.cond = dfsin.flows ./ dfsin.toflux .* 100
dfsin = sort(dfsin, [:todist, :cond], rev = true)
dfsin.id = repeat(1:399, 400)

function condsum(df, type)
    col = type == "outflux" ? :fromdist : :todist
    Norigins = length(unique(df[!, col]))
    if Norigins == 1
        dfss = combine(groupby(df, :id), :cond => mean,
                       col => first => col)
    else
        dfss = combine(groupby(df, :id), :cond => mean)
    end
    dfss.cond_sum = cumsum(dfss.cond_mean)
    return dfss
end

condout = condsum(dfsout, "outflux")
condin = condsum(dfsin, "influx")

f = Figure(size = (400, 400), fontsize = 10);
ax1 = Axis(f[1, 1:2];
           title = "How much does each percentile\ncontribute to total moves?",
           subtitle = "2000-2017, all age groups",
           xlabel = "flow percentile\n(flow value)",
           ylabel = "% of total moves",
           xticks = (ticks, vs),
           xgridvisible = false, ygridvisible = false)
barplot!(ax1, df1.qs, df1.flows)

ax2 = Axis(f[2, 1],
           title = "How many destinations are needed\nto capture x% of moves?",
           subtitle = "2000-2017, all age groups, all origins",
           xlabel = "Top Destinations (top to bottom)",
           ylabel = "Cumulated Relative Flows",
           xgridvisible = false, ygridvisible = false)
Makie.xlims!(ax2, (0, 100))
Makie.ylims!(ax2, (0, 100))

for o in unique(dfsout.fromdist)
    foo = condsum(origin(dfsout, o), "outflux")
    lines!(ax2, foo.id, foo.cond_sum, color = :lightgrey, alpha = .2)
end
dfss = condsum(dfsout, "outflux")
lines!(ax2, dfss.id, dfss.cond_sum)

ax3 = Axis(f[2, 2],
           title = "How many origins are needed\nto capture x% of moves?",
           subtitle = "2000-2017, all age groups, all destinations",
           xlabel = "Top Origins (top to bottom)",
           ylabel = "Cumulated Relative Flows",
           xgridvisible = false, ygridvisible = false)
Makie.xlims!(ax3, (0, 100))
Makie.ylims!(ax3, (0, 100))

for o in unique(dfsin.todist)
    foo = condsum(destination(dfsin, o), "influx")
    lines!(ax3, foo.id, foo.cond_sum, color = :lightgrey, alpha = .2)
end
dfss = condsum(dfsin, "influx")
lines!(ax3, dfss.id, dfss.cond_sum)
save(joinpath(outp, "quantiles_top.pdf"), f)

foo = reduce(vcat, [condsum(origin(dfs, o)) for o in unique(dfs.fromdist)])
leftjoin!(foo, diy[!, [:distcode, :name]], on = :fromdist => :distcode)
foo = sort(foo[foo.id .== 50, :], :fromdist)

fig = Figure(size = (400, 400), fontsize = 10);
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title = "Flow Dispersion",
          subtitle = "All age group, yearly average 2000 - 2017")
viz!(ax, rmreg(shp, :AGS).geometry, color = foo.cond_sum, colorrange = extrema(foo.cond_sum));
overlay_states(ax, st)
Colorbar(fig[2, 1], limits = extrema(foo.cond_sum), vertical = false,
         height = 3, width = Relative(.4), label = "Dispersion")
fig
