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
out = combine(groupby(df, [:fromdist, :year, :agegroup]), :flows => sum => :outflow)
leftjoin!(df, out, on = [:fromdist, :year, :agegroup])
df.cond = df.flows ./ df.outflow .* 100

dfs = combine(groupby(df, [:fromdist, :todist]),
              [:flows => sum => :flows])
leftjoin!(dfs, outflux(df, sum), on = :fromdist)
dfs.cond = dfs.flows ./ dfs.outflux .* 100
dfs = sort(dfs, [:fromdist, :cond], rev = true)
dfs.id = repeat(1:399, 400)

function condsum(df)
    Norigins = length(unique(df.fromdist))
    if Norigins == 1
        dfss = combine(groupby(df, :id), :cond => mean,
                       :fromdist => first => :fromdist)
    else
        dfss = combine(groupby(df, :id), :cond => mean)
    end
    dfss.cond_sum = cumsum(dfss.cond_mean)
    return dfss
end

f = Figure(size = (400, 300), fontsize = 10);
ax1 = Axis(f[1, 1];
           title = "How much does each percentile\ncontribute to total moves?",
           subtitle = "2000-2017, all age groups",
           xlabel = "flow percentile\n(flow value)",
           ylabel = "% of total moves",
           xticks = (ticks, vs),
           xgridvisible = false, ygridvisible = false)
barplot!(ax1, df1.qs, df1.flows)

ax2 = Axis(f[1, 2],
           title = "How many destinations are needed\nto capture x% of moves?",
           subtitle = "2000-2017, all age groups, all origins",
           xlabel = "Top Destinations (top to bottom)",
           ylabel = "Cumulated Relative Flows",
           xgridvisible = false, ygridvisible = false)
Makie.xlims!(ax2, (0, 100))
Makie.ylims!(ax2, (0, 100))

for o in unique(dfs.fromdist)
    foo = condsum(origin(dfs, o))
    lines!(ax2, foo.id, foo.cond_sum, color = :lightgrey, alpha = .2)
end
dfss = condsum(dfs)
lines!(ax2, dfss.id, dfss.cond_sum)
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
