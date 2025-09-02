using DataFrames, CSV, CairoMakie, Statistics, StatsBase, Revise
using LaTeXStrings, AlgebraOfGraphics
outp = "/home/konstantin/paper/sections/texdocs/images/"

includet("../src/diagplots.jl")

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
years = unique(df.year)
ages = unique(df.agegroup)


######################################################################
################# total and relative flows ###########################
dft = combine(groupby(df, [:year, :agegroup]), :flows => sum => :total)
leftjoin!(dft, pop, on = [:year, :agegroup])
dft.prob = dft.total ./ dft.agepop .* 100
dft.total = dft.total ./ 100e3

f = Figure();
ax1 = Axis(f[1, 1], title = "Total Flows", ylabel = "Flows (100k)", xlabel = "Year")
ax2 = Axis(f[1, 2], ylabel = "Movement Frequency (%)", xlabel = "Year")
p1 = draw!(ax1, data(dft) * mapping(:year, :total, color = :agegroup) * visual(Lines))
p2 = draw!(ax2, data(dft) * mapping(:year, :prob, color = :agegroup) * visual(Lines))
legend!(f[1, 3], p1)
save(joinpath(outp, "flows.pdf"), f)

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
    df1 = reduce(vcat, [condsum(df, y, a) for y in years])
    df1 = combine(groupby(df1, [:agegroup, :id]), :sum => mean => :sum)
    return df1
end
df1 = reduce(vcat, [condsumage(df, a) for a in ages])
p1 = plot(df1.id, df1.sum, group = df1.agegroup,
          xlim = (0, 100), ylim = (0, 100),
          xlab = "Most important destinations",
          ylab = "cumulated relative flows",
          title = "How many destinations are needed\nto capture x % of flows?")
savefig(p1, joinpath(outp, "topdest.pdf"))


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
####################### Flow Magnitude ###############################
df1 = sort(outflux(df, sum, [:todist]), :outflux)
flows = df1.outflux ./ length(years)

qs = quantile(flows, range(0, 1, 100))
flows_sum = [sum(flows[flows .> qs[i-1] .&& flows .< qs[i]]) for i in 2:length(qs)]
ticks = [.25, .5, .75, .9, .95, 1]
vs = Int.(round.(quantile(flows, ticks), digits = 0))
vs = string.(Int.(ticks .* 100)) .* "\n" .* "(" .* string.(vs) .* ")"

p1 = bar(flows_sum ./ sum(flows) .* 100, ylim = (0, 50), label = "",
         xlab = "flow percentile, (flow value)", ylab = "% of total flows",
         title = "How much does each percentile\ncontribute to total flows?",
         xticks = (ticks .* 100, tick_vals))
savefig(p1, joinpath(outp, "flowquantiles.pdf"))

quantile(flows, .5)
1 - sum(flows[flows .<= quantile(flows, .99)]) ./ sum(flows)
1 - sum(flows[flows .<= quantile(flows, .9)]) ./ sum(flows)
sum(flows[flows .<= quantile(flows, .75)]) ./ sum(flows)
