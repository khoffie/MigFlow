using DataFrames, CSV, CairoMakie, Statistics, StatsBase, Revise
using LaTeXStrings, AlgebraOfGraphics, GeoIO, GeoStats
outp = "/home/konstantin/paper/sections/texdocs/images/"

includet("../src/diagplots.jl")
includet("../src/analyzegeo.jl")
includet("../src/model_helpers.jl")

age(df, age) = filter(:agegroup => n -> n ∈ age, df)
year(df, y) = filter(:year => n -> n ∈ y, df)
origin(df, o) = filter(:fromdist => n -> n ∈ o, df)
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
sort(net[!, [:yearly_total, :name, :fromdist]], :yearly_total)

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

function plot()
    f = Figure(size = (300, 300), fontsize = 10);
    ax1 = Axis(f[1, 1], title = "Total Flows", ylabel = "Flows (100k)", xlabel = "Year");
    ax2 = Axis(f[1, 2], title = "Relative Frequency", ylabel = "Movement Frequency (%)", xlabel = "Year");
    p1 = draw!(ax1, data(dft) * mapping(:year, :total, color = :agegroup) * visual(Lines));
    p2 = draw!(ax2, data(dft) * mapping(:year, :prob, color = :agegroup) * visual(Lines));
    legend!(f[0, 1:2], p1; framevisible = true, rowgap = 0, colgap = 0, orientation = :horizontal,
            markersize = 5, padding = (0, 0, 0, 0), patchlabelgab = 10, patchsize = (10, 10))
    return f
end

save(joinpath(outp, "flows.pdf"), plot())

######################################################################
####################### Flow Magnitude ###############################
df1 = sort(outflux(df, sum, [:todist]), :outflux)
flows = df1.outflux ./ length(years)

qs = quantile(flows, range(0, 1, 100))
flows_sum = [sum(flows[flows .> qs[i-1] .&& flows .< qs[i]]) for i in 2:length(qs)]
flows_sum = flows_sum ./ sum(flows) .* 100
df1 = DataFrame(qs = 1:(length(qs) - 1), flows = flows_sum)
ticks = [.25, .5, .75, .9, .95, 1] * 100
vs = Int.(round.(quantile(flows, ticks ./ 100), digits = 0))
vs = string.(Int.(ticks)) .* "\n" .* "(" .* string.(vs) .* ")"

1 - sum(flows[flows .<= quantile(flows, 0.99)]) / sum(flows)
1 - sum(flows[flows .<= quantile(flows, 0.9)]) / sum(flows)
sum(flows[flows .<= quantile(flows, 0.5)]) / sum(flows)


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
