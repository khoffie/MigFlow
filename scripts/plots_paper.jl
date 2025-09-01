using DataFrames, CSV, StatsPlots
outp = "/home/konstantin/paper/sections/texdocs/images/"

df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
di = CSV.read("../data/districts.csv", DataFrame)
pop = combine(groupby(unique(df, [:fromdist, :year, :agegroup]),
                      [:year, :agegroup]), :frompop => sum => :agepop)
## combine(groupby(pop, :year), :agepop => sum)


## total and relative flows
dft = combine(groupby(df, [:year, :agegroup]), :flows => sum => :total)
leftjoin!(dft, pop, on = [:year, :agegroup])
dft.prob = dft.total ./ dft.agepop .* 100

p1 = plot(dft.year, dft.total, group = dft.agegroup, linewidth = 2,
          ylab = "total flows", label = "")
p2 = plot(dft.year, dft.prob, group = dft.agegroup, linewidth = 2,
          ylab = "probability to move")
savefig(plot(p1, p2), joinpath(outp, "flows.pdf"))
