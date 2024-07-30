using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random, ReverseDiff, Revise, RCall
using OptimizationOptimJL
includet("model2.jl")
includet("Rutils.jl")
includet("simulateddata.jl")

Random.seed!(20240719)

## read in data, it should have the following columns:
## flows, fromdist, todist, frompop, topop, distance, gdpcfrom, gdpcto, agegroup

sims = false
if sims
#    ourdat = CSV.read("data/simulations.csv",DataFrame)
#    rename!(ourdat,Dict("predict" => "flows"))
#    ourdat.flows = round.(Int32,ourdat.flows)
    ourdat,knowndesir = simdatafromtemplate(Xoshiro(20240725),"data/flowtemplate.csv")
    rename!(ourdat,Dict(:age_group => :agegroup,:distance => :dist,:flow => :flows))
else
    ourdat = CSV.read("/home/donkon/Diss/inst/extdata/clean/daniel/FlowData.csv",DataFrame)
end
ourdat.fromdist = categorical(ourdat.fromdist)
ourdat.todist = categorical(ourdat.todist)
ourdat.agegroup = categorical(ourdat.agegroup)
levels!(ourdat.agegroup,["below18","18-25","25-30","30-50","50-65","above65"])
rename!(ourdat, Dict(:dist => :distance))

meddist = median(ourdat.distance)
Ndist = length(unique(ourdat.fromdist))
Nages = length(unique(ourdat.agegroup))
##length(unique(levelcode.(ourdat.agegroup)))

ourdat2 = ourdat
#ourdat2 = ourdat[sample(1:nrow(ourdat),20000),:]
## inis = 1 .+ 0.05 .* randn(Ndist)
ourdat2

netactual = calcnet(ourdat2.flows,
                    levelcode.(ourdat2.fromdist),
                    levelcode.(ourdat2.todist),
                    levelcode.(ourdat2.agegroup),Nages,Ndist)

model2 = migration2(ourdat2.flows,levelcode.(ourdat2.fromdist),levelcode.(ourdat2.todist),ourdat2.frompop,ourdat2.topop,ourdat2.distance, levelcode.(ourdat2.agegroup), Nages, Ndist, meddist,netactual)

opinit = [[11.0,3.3,1.8,.1,.5] ; ones(Ndist*6)];
opinit = [[10.40,3.1,1.25,.051,.25] ; repeat(knowndesir,6)];

## use optimization to find a good fit
mapfit2 = maximum_a_posteriori(model2, LBFGS() ; adtype = AutoReverseDiff(), 
            initial_params = opinit, maxiters = 200, maxtime = 600, reltol = .08,
            lb = zeros(length(opinit)), ub=25*ones(length(opinit)))
mapfit2
initvals = mapfit2.values

df = DataFrame(param = names(initvals, 1), estim = values(initvals))
CSV.write("./data/opti_vals.csv", df)

plotdesirability(initvals)
plotnetmigration(netmigr)

param_values = values(mapfit2.values)

println("Optimal values are: ")
@show initvals

initvals2 = Iterators.Repeated(initvals .+ rand(Uniform(0.0,.50),length(initvals)))

initvals2 = Iterators.Repeated(opinit)

println("Perturbed Initial values are:")

@show iterate(initvals2)[1]
## start the sampling at a location biased away from the mode, by increasing all parameters 
## by a small uniform perturbation (this avoids anything that has to be positive becoming negative)

fit2 = sample(model2, NUTS(100,.8; adtype=AutoReverseDiff(true)),
              MCMCThreads(), 100, 3; 
              init_params = initvals2)

display(fit2)

mainparms2 = fit2[:,[:a,:b,:c,:d0,:neterr],:]
plot(mainparms2) |> display()

