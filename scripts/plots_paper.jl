using DataFrames, CSV, CairoMakie, Statistics, StatsBase, Revise
using LaTeXStrings, AlgebraOfGraphics, GeoIO, GeoStats, GLM
using KernelDensity
outp = "/home/konstantin/paper/sections/texdocs/images/"

includet("../src/diagplots.jl")
includet("../src/analyzegeo.jl")
includet("../src/model_helpers.jl")

age(df, age::Vector{String}) = filter(:agegroup => n -> n ∈ age, df)
age(df, age::String) = filter(:agegroup => n -> n == age, df)
year(df, y::Vector{Int64}) = filter(:year => n -> n ∈ y, df)
year(df, y::Int64) = filter(:year => n -> n == y, df)
origin(df, o::Vector{Int64}) = filter(:fromdist => n -> n ∈ o, df)
origin(df, o::Int64) = filter(:fromdist => n -> n == o, df)
destination(df, d) = filter(:todist => n -> n ∈ d, df)
code(df, c) = filter(:distcode => n -> n ∈ c, df)

function topn(df, group, col, N = 10)
    dfg = groupby(sort(df, col, rev = true), group)
    dfg = [dfg[i][1:N, [group..., col]] for i in 1 : length(dfg)]
    return reduce(vcat, dfg)
end

function outflux(df, f, group = nothing)
    group == nothing ? g = :fromdist : g = [:fromdist, group...]
    out = combine(DataFrames.groupby(df, g), :flows => f)
    return rename!(out, "flows" * "_$(string(f))" => :outflux)
end

df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
di = CSV.read("../data/districts.csv", DataFrame)
pop = combine(groupby(unique(df, [:fromdist, :year, :agegroup]),
                      [:year, :agegroup]), :frompop => sum => :agepop)
## combine(groupby(pop, :year), :agepop => sum)
shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp")

years = unique(df.year)
ages = unique(df.agegroup)


######################################################################
################# total flows per region  ############################
net = calc_net(df, :flows);
net.geo = log.(net.total ./ length(years));
net.yearly_total = net.total ./ length(years)
leftjoin!(net, year(di, 2017)[!, [:distcode, :name]], on = :fromdist => :distcode)

fig = Figure(size = (400, 400), fontsize = 10)
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title = "Which regions contribute the most to moves?",
          subtitle = "All age group, yearly average 2000 - 2017")
tightlimits!(ax)
hidedecorations!(ax)
hidespines!(ax)
viz!(ax, shp.geometry, color = net.geo, colorrange = extrema(net.geo));
overlay_states(ax, st)
Colorbar(fig[2, 1], limits = extrema(net.geo), vertical = false,
         height = 3, width = Relative(.4), label = "log(influx + outflux)")
resize_to_layout!(fig)
save(joinpath(outp, "totalmap.pdf"), fig)
sum(origin(net, [:2000, :9162, :11000]).yearly_total) / 2.2e6

######################################################################
################# total and relative flows ###########################
dft = combine(groupby(df, [:year, :agegroup]), :flows => sum => :total)
leftjoin!(dft, pop, on = [:year, :agegroup])
dft.prob = dft.total ./ dft.agepop .* 100
dft.total = dft.total ./ 100e3

f = Figure(size = (300, 300), fontsize = 10);
ax1 = Axis(f[1, 1],
           ylabel = "Flows (100k)",
           xlabel = "Year",
           xgridvisible = false, ygridvisible = false);
ax2 = Axis(f[1, 2],
           ylabel = "Movement Frequency (%)",
           xlabel = "Year",
           xgridvisible = false, ygridvisible = false);
for a in ages
    df = age(dft, [a])
    lines!(ax1, df.year, df.total, label = a)
    lines!(ax2, df.year, df.prob, label = a)

end
axislegend(ax2, framevisible = false, patchsize = (10, 10), position = (.3, .5))
save(joinpath(outp, "flows_age.pdf"), f)

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

function condsum(df, y, a)
    dfe = sort(age(year(df, y), a)[
        !, [:year, :agegroup, :fromdist, :todist, :cond]],
               [:fromdist, :cond], rev = true)
    dfe.id = repeat(1:400, 401)
    df2 = combine(groupby(dfe, :id), :cond => mean)
    df2.agegroup .= unique(dfe.agegroup)[1]
    df2.year .= unique(dfe.year)[1]
    df2.sum = cumsum(df2.cond_mean)
    return df2
end

function condsumage(df, a)
    df1 = reduce(vcat, [condsum(df, y, [a]) for y in years])
    df1 = combine(groupby(df1, [:agegroup, :id]), :sum => mean => :sum)
    return df1
end
df2 = reduce(vcat, [condsumage(df, a) for a in ages])
df2 = combine(groupby(df2, :id), :sum => mean => "sum")

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
xlims!(ax2, (0, 100))
ylims!(ax2, (0, 100))
lines!(ax2, df2.id, df2.sum)
save(joinpath(outp, "quantiles_top.pdf"), f)

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

######################################################################
##################### Data Issues ####################################
df1 = destination(origin(df, [3159]), [5978])
combine(groupby(df1, :year), :flows => sum)

sort(outflux(origin(df, 3159), [:year]), :outflux)
sort(outflux(destination(origin(df, 3159), 5978), [:year]), :outflux)

######################################################################
############################# GLM ####################################

df17 = combine(groupby(year(df, 2017), [:fromdist, :todist]),
               [:flows => sum => :flows,
                :topop => first => :topop,
                :dist => first => :dist,
               :frompop => sum => :frompop])

sum(unique(df17, :topop).topop)
sum(unique(df17, :frompop).frompop)

m = glm(@formula(flows ~ log(dist)), df17, Poisson(), LogLink())
m = glm(@formula(flows ~ log(frompop) + log(topop) + log(dist)), df17, Poisson(), LogLink())
1 - var(df17.flows .- predict(m)) / var(df17.flows)

######################################################################
########################### Distance #################################
res, dfres = maincircle(df17, di, false)
res
dfres = filter(:fpp => !isnan, dfres)
dfres.distbin = floor.(dfres.dist ./ 20) .* 20  # Group into bins like [0,1), [1,2), etc.
dfres.fppop = dfres.flows ./ dfres.topop .* 100000
leftjoin!(dfres, year(di, 2017)[!, [:distcode, :area]], on = :todist => :distcode)
dfres.fpp = dfres.flows ./ (dfres.frompop ./ dfres.area )


dfres2 = combine(groupby(dfres, :distbin), [:flows, :fpp, :fppop] .=> sum)
y1 = log.(dfres2.flows_sum)
y2 = log.(dfres2.fpp_sum)
y3 = log.(dfres2.fppop_sum)

f = Figure(size = (350, 350), fontsize = 10);
ax = Axis(f[1, 1],
          xlabel = "Distance (km)",
          xgridvisible = false, ygridvisible = false)
hideydecorations!(ax)
xs = StatsBase.wsample(dfres.dist, Weights(dfres.pot), 10^4)
lines!(ax, kde(xs, Gamma(1, 10)), label = L"\text{Potential} ∝ a2πr", color = :purple)
xs = StatsBase.sample(df.dist, 10^4)
lines!(ax, kde(xs, Gamma(1, 10)), label = L"\text{Potential} ∝ Pop", color = :green)
xs = StatsBase.wsample(df.dist, Weights(df.flows), 10^4)
lines!(ax, kde(xs, Gamma(1, 10)), label = L"\text{Observed}", color = :red)
axislegend(ax, framevisible = false, patchsize = (5, 10))

ax2 = Axis(f[1, 2],
           xlabel = "Distance (km)",
           xgridvisible = false, ygridvisible = false)
hideydecorations!(ax2)

xs = StatsBase.wsample(dfres.dist, Weights(dfres.fpp), 10^4)
lines!(ax2, kde(xs, Gamma(1, 10)), label = L"\text{Adjusted} a2πr", color = :purple)
xs = StatsBase.wsample(dfres.dist, Weights(dfres.fppop), 10^4)
lines!(ax2, kde(xs, Gamma(1, 10)), label = L"\text{Adjusted} Pop", color = :green)
xs = StatsBase.wsample(dfres.dist, Weights(dfres.flows), 10^4)
lines!(ax2, kde(xs, Gamma(1, 10)), label = L"\text{Observed}", color = :red)
axislegend(ax2, framevisible = false, patchsize = (5, 10))

ax3 = Axis(f[1, 3],
           xlabel = "Distance (km)",
           ylabel = "Density of log(y)",
           xgridvisible = false, ygridvisible = false)
hideydecorations!(ax3)
lines!(ax3, dfres2.distbin, y2 ./ sum(y2), color = :purple,
       label = "Adjusted\ny = flows / a2πr"
       )
lines!(ax3, dfres2.distbin, y3 ./ sum(abs.(y3)), color = :green,
       label = "Adjusted\ny = flows / P")
lines!(ax3, dfres2.distbin, y1 ./ sum(y1), color = :red,
       label = "Observed\ny = flows")
axislegend(ax3, framevisible = false, patchsize = (5, 10))

titlelayout = GridLayout(f[0, 1], halign = :left, tellwidth = false)
Label(titlelayout[1, 1], "Movement Distances", halign = :left, fontsize = 15, font = "TeX Gyre Heros Bold Makie")
Label(titlelayout[2, 1], "2000-2017, all age groups", halign = :left, fontsize = 10)
rowgap!(titlelayout, 0)
f
save(joinpath(outp, "distances.pdf"), f)
StatsPlots.histogram(year(di, 2017).density)

sum(destination(dfm, 11000).outflux) / sum(dfm.outflux)

test = origin(df, 9663)

destination(origin(df, 9663), 11000)

t = test[test.dist .> 300 .&& test.dist .< 400, :]
t = outflux(t, sum, [:todist])

sum(destination(t, 11000).flows) / sum(t.flows)
leftjoin!(t, year(di, 2017)[!, [:distcode, :area, :pop]], on = :todist => :distcode)

t.y = t.outflux ./ t.area
t
destination(t, 11000)

mean(log.(t.y))

Plots.scatter(log.(t.pop), log.(t.outflux))
t
names(di)
di.area = di.pop ./ di.density
StatsPlots.scatter(log.(year(di, 2017).area), log.(year(di, 2017).pop))

code(year(di, 2017), 11000)
exp(7)
log(900)
df[df.dist .> 800, :]
