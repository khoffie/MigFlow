using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random, ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions, ApproxFun, Serialization, Printf, DataFramesMeta
## includet("Rutils.jl")
## includet("simulateddata.jl")
includet("debughelpers.jl")

Random.seed!(20240719)

munis = CSV.read("./data/munis_pop.csv", DataFrame)
rename!(munis,Dict(:age_group => :agegroup))
coords = CSV.read("./data/district_coords.csv", DataFrame)
sims = true
if sims
#    dt = CSV.read("data/simulations.csv",DataFrame)
#    rename!(dt,Dict("predict" => "flows"))
#    dt.flows = round.(Int32,dt.flows)
    dt,knowndesir = simdatafromtemplate(Xoshiro(20240725),"data/flowtemplate.csv")
    rename!(dt,Dict(:age_group => :agegroup,:distance => :dist,:flow => :flows))
else
    dt = CSV.read("/home/konstantin/Diss/inst/extdata/clean/daniel/FlowData.csv",DataFrame)
end
dt2 = dt
# dt2 = dt[dt[: , 6] .> 10, :]
# random40 = selectdists(munis, 50)[11:50]
# top40 = selectdists(munis, 400)[1:100]
# dt2 = subsetdists(dt, top40)

## top10 from R: 11000  2000  9162  3241  5315  6412  8111  5111  5562  5382
dt2.fromdist = categorical(dt2.fromdist)
dt2.todist = categorical(dt2.todist)
dt2.agegroup = categorical(dt2.agegroup)
levels!(dt2.agegroup,["below18","18-25","25-30","30-50","50-65","above65"])
rename!(dt2, Dict(:dist => :distance))
## scatter(dt2.distance, log.(dt2.flows ./ dt2.frompop ./ dt2.topop))

meddist = median(dt2.distance)
Ndist = length(unique(dt2.fromdist))
Nages = length(unique(dt2.agegroup))

