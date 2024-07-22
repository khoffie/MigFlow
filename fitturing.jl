using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random, ReverseDiff
include("model2.jl")
Random.seed!(20240719)

## read in data, it should have the following columns:
## flows, fromdist, todist, frompop, topop, distance, gdpcfrom, gdpcto, agegroup
ourdat = CSV.read("/home/donkon/Diss/inst/extdata/clean/daniel/FlowData.csv",DataFrame)
ourdat.fromdist = categorical(ourdat.fromdist)
ourdat.todist = categorical(ourdat.todist)
ourdat.agegroup = categorical(ourdat.agegroup)
levels!(ourdat.agegroup,["below18","18-25","25-30","30-50","50-65","above65"])
rename!(ourdat, Dict(:dist => :distance))

meddist = median(ourdat.distance)
meddist = 293.0 ## median distance between districts according to Konstantin's printed report in km
Ndist = length(unique(ourdat.fromdist))
Nages = length(unique(ourdat.agegroup))
##length(unique(levelcode.(ourdat.agegroup)))


ourdat2 = ourdat[sample(1:nrow(ourdat),1000),:]
inis = 1 .+ 0.05 .* randn(Ndist)

fit2 = sample(migration2(ourdat2.flows,levelcode.(ourdat2.fromdist),levelcode.(ourdat2.todist),ourdat2.frompop,ourdat2.topop,ourdat2.distance, levelcode.(ourdat2.agegroup), Nages, Ndist, meddist,netactual),
              NUTS(500,.8; adtype=AutoReverseDiff(true)),
              MCMCThreads(), 600, 3)

# ; init_params = Iterators.repeated(inis)
## too many parameters for this plot
#plot(fit2) |> display()

display(fit2)

mainparms2 = fit2[:,[:a,:b,:c,:d0,:neterr],:]
plot(mainparms2) |> display()
