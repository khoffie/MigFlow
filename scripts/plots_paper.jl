using DataFrames, CSV, CairoMakie, Statistics, StatsBase, Revise
using LaTeXStrings, AlgebraOfGraphics, GeoIO, GeoStats, GLM
using KernelDensity
outp = "/home/konstantin/paper/sections/texdocs/images/"

includet("../src/diagplots.jl")
includet("../src/analyzegeo.jl")
includet("../src/model_helpers.jl")
includet("../src/samplecircle.jl")

age(df, age::Vector{String}) = filter(:agegroup => n -> n ∈ age, df)
age(df, age::String) = filter(:agegroup => n -> n == age, df)
year(df, y::Vector{Int64}) = filter(:year => n -> n ∈ y, df)
year(df, y::Int64) = filter(:year => n -> n == y, df)
origin(df, o::Vector{Int64}) = filter(:fromdist => n -> n ∈ o, df)
origin(df, o::Int64) = filter(:fromdist => n -> n == o, df)
destination(df, d) = filter(:todist => n -> n ∈ d, df)
code(df, c) = filter(:distcode => n -> n ∈ c, df)
rmreg(df::DataFrame, col) = filter(col => n -> n != 3159, df) ## data issues
rmreg(df::GeoTable, col) = GeoTable(filter(col => n -> n != 3159, DataFrame(df))) ## data issues

function hideall!(ax)
    tightlimits!(ax)
    hidedecorations!(ax)
    hidespines!(ax)
end

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
di.area = di.pop ./ di.density;

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
net = rmreg(net, :fromdist)
net.geo = log.(net.total ./ length(years));
net.yearly_total = net.total ./ length(years)

diy = combine(groupby(di, :distcode), [
    ## should pop be pop of Germans only?
    :pop => mean => :pop,
    :name => first => :name
])

leftjoin!(net, diy, on = :fromdist => :distcode)
net.flow_share = net.yearly_total ./ net.pop .* 100

fig = Figure(size = (400, 400), fontsize = 10);
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title = "Influx + Outflux",
          subtitle = "All age group, yearly average 2000 - 2017")
viz!(ax, rmreg(shp, :AGS).geometry, color = net.geo, colorrange = extrema(net.geo));
overlay_states(ax, st)

ax2 = Axis(fig[1, 2], aspect=DataAspect(),
          title = "(Influx + Outflux) / Population",
          subtitle = "All age group, yearly average 2000 - 2017")
viz!(ax2, rmreg(shp, :AGS).geometry, color = net.flow_share, colorrange = extrema(net.flow_share));
overlay_states(ax2, st)

Colorbar(fig[2, 1], limits = extrema(net.geo), vertical = false,
         height = 3, width = Relative(.4), label = "log(influx + outflux)")
Colorbar(fig[2, 2], limits = extrema(net.flow_share), vertical = false,
         height = 3, width = Relative(.4), label = "(influx + outflux) / pop")

hideall!(ax)
hideall!(ax2)
fig
resize_to_layout!(fig)
save(joinpath(outp, "totalmap.pdf"), fig)




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
leftjoin!(df17, year(di, 2017)[!, [:distcode, :area]],
          on = :todist => :distcode)
df17.toarea = Float64.(df17.area)
select!(df17, Not(:area))
leftjoin!(df17, year(di, 2017)[!, [:distcode, :area]],
          on = :fromdist => :distcode)
df17.fromarea = Float64.(df17.area)
select!(df17, Not(:area))

sum(unique(df17, :topop).topop)
sum(unique(df17, :frompop).frompop)

m = glm(@formula(flows ~ log(dist)), df17, Poisson(), LogLink())
m = glm(@formula(flows ~ log(frompop) + log(topop) + log(dist)), df17, Poisson(), LogLink())
1 - var(df17.flows .- predict(m)) / var(df17.flows)

######################################################################
########################### Distance #################################
res = maincircle(df17, di, false);
dfdist = combine(groupby(df17, :dist),
                 [ nrow => :count,
                   :flows => sum => :flows,
                   [:topop, :frompop] => ((x, y) -> sum(log.(x) + log.(y))) => :pop,
                   [:toarea, :fromarea] => ((x, y) -> sum(log.(x) + log.(y))) => :area
                   ])
res = dropmissing(leftjoin(res, dfdist, on = :radius => :dist))

model(radius, γ, ϵ) = Float64.(radius).^γ .+ ϵ

function plot(res, γ, ϵ)
    res.pred = model(res, γ, ϵ)

    f = Figure(size = (350, 350), fontsize = 10);
    ax1 = Axis(f[1, 1],
              xlabel = "Distance (km)",
              xgridvisible = false, ygridvisible = false)

    xs = StatsBase.sample(res.radius, Weights(res.pop), 10^4);
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :red, label = "Pop")
    xs = StatsBase.sample(res.radius, Weights(res.area), 10^4);
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :blue, label = "Area")
    xs = StatsBase.sample(res.radius, Weights(res.pot), 10^4)
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :purple, label = "a2πr")
    xs = StatsBase.sample(res.radius, Weights(res.count), 10^4)
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :green, label = "Unweighted")
    xs = StatsBase.sample(res.radius, Weights(res.flows), 10^4)
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :yellow, label = "Flows")
    xs = StatsBase.sample(res.radius, Weights(res.count .* res.pred), 10^4)
    lines!(ax1, kde(xs, Gamma(1, 10)), color = :brown, label = "radius^$γ + $ϵ")

    ax2 = Axis(f[1, 2],
               xlabel = "Distance (km)",
               xgridvisible = false, ygridvisible = false)
    res.flows[res.flows .== 0] .= 1

    xs = StatsBase.sample(res.radius, Weights(res.flows ./ res.flows), 10^4);
    lines!(ax2, kde(xs, Gamma(1, 10)), color = :red, label = "Flows / Flows")
    xs = StatsBase.sample(res.radius, Weights(res.flows ./ (res.pred .* res.count)), 10^4);
    lines!(ax2, kde(xs, Gamma(1, 10)), color = :brown, label = "Flows / (radius^$γ + $ϵ)")

    hideydecorations!(ax1)
    hideydecorations!(ax2)
    axislegend(ax1, framevisible = false, patchsize = (5, 10))
    axislegend(ax2, framevisible = false, patchsize = (5, 10))

    titlelayout = GridLayout(f[0, 1], halign = :left, tellwidth = false)
    Label(titlelayout[1, 1], "Movement Distances", halign = :left, fontsize = 15, font = "TeX Gyre Heros Bold Makie")
    Label(titlelayout[2, 1], "2000-2017, all age groups", halign = :left, fontsize = 10)
    rowgap!(titlelayout, 0)

    return res, f
end

r, f = plot(res, (-2.05), 10^(-5))
save(joinpath(outp, "distances.pdf"), f)
