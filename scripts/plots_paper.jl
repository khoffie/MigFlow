using DataFrames, CSV, StatsPlots, Statistics, StatsBase
outp = "/home/konstantin/paper/sections/texdocs/images/"

age(df, age) = filter(:agegroup => n -> n == age, df)
year(df, y) = filter(:year => n -> n == y, df)
origin(df, o) = filter(:fromdist => n -> n ∈ o, df)
function topn(df, group, col, N = 10)
    dfg = groupby(sort(df, col, rev = true), group)
    dfg = [dfg[i][1:N, [group..., col]] for i in 1 : length(dfg)]
    return reduce(vcat, dfg)
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

p1 = plot(dft.year, dft.total, group = dft.agegroup, linewidth = 2,
          ylab = "total flows", label = "")
p2 = plot(dft.year, dft.prob, group = dft.agegroup, linewidth = 2,
          ylab = "probability to move")
savefig(plot(p1, p2), joinpath(outp, "flows.pdf"))


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
