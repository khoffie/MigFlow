using CSV, DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random, ReverseDiff, Revise, RCall
using OptimizationOptimJL, Distributions
includet("model2.jl")
includet("Rutils.jl")
includet("simulateddata.jl")
includet("debughelpers.jl")

Random.seed!(20240719)

munis = CSV.read("./data/munis_pop.csv", DataFrame)
rename!(munis,Dict(:age_group => :agegroup))
sims = true
if sims
#    dt = CSV.read("data/simulations.csv",DataFrame)
#    rename!(dt,Dict("predict" => "flows"))
#    dt.flows = round.(Int32,dt.flows)
    dt,knowndesir = simdatafromtemplate(Xoshiro(20240725),"data/flowtemplate.csv")
    rename!(dt,Dict(:age_group => :agegroup,:distance => :dist,:flow => :flows))
else
    dt = CSV.read("./data/FlowDataGermans.csv",DataFrame)
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

netactual = calcnet(dt2.flows,
                    levelcode.(dt2.fromdist),
                    levelcode.(dt2.todist),
                    levelcode.(dt2.agegroup),Nages,Ndist)

model2 = migration2(dt2.flows, levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                    dt2.frompop, dt2.topop, dt2.distance,
                    levelcode.(dt2.agegroup), Nages, Ndist, meddist, netactual)


model2
opinit = [fill(11.0,Nages); fill(3.3,Nages); fill(1.8,Nages); fill(.1,Nages); [1]; 1 .* 100 .* ones(Ndist*Nages)]

opinit = [  [.5, 4.5, 5.7, .53, .2, .14];
            [5.6, 7.1, 7.9, 5.8, 3.8, 3.7];
            [4.0, 4.0, 4.0, 4.5, 4.0, 4.0]; 
            [.6, .8, .75, .64, .5, .47]; 
            [1.0]; 
            1 .* 100 .* ones(Ndist*Nages)]

lower = [fill(0,Nages); fill(0,Nages); fill(0,Nages); fill(0,Nages); [0]; .7 .* 100 .* ones(Ndist*Nages)]
upper = [fill(20,Nages); fill(20,Nages); fill(5,Nages); fill(1,Nages); [2]; 1.3 .* 100 .* ones(Ndist*Nages)]

# appx = vi(model2, ADVI(10, 1000))
# rand(appx, 1000)
## use optimization to find a good fit
##lb = lower, ub = upper,
mapfit2 = maximum_a_posteriori(model2, LBFGS() ; adtype = AutoReverseDiff(), 
                               initial_params = opinit, lb = lower, ub = upper,
                               maxiters = 200, maxtime = 60, reltol = .008)


## started 14:59
mapfit2
initvals = mapfit2.values
initvals2 = Iterators.Repeated(initvals .+ rand(Uniform(0.0,.50),length(initvals)))
# df = DataFrame(param = names(initvals, 1), estim = values(initvals))
# CSV.write("./data/opti_vals.csv", df)

fit2 = sample(model2, NUTS(500,.8; adtype=AutoReverseDiff(true)), 100,
              init_params = opinit,
              verbose = true, progress = true)
fit2
plot(fit2)
display(fit2)

mainparms2 = fit2[:,[:a,:b,:c,:d0,:neterr],:]
plot(mainparms2) |> display()


opinit

chain = sample(model2, Prior(), 30)
chain
plot(chain)
loglikelihood(model2, opinit)
