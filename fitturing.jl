
using CSV,DataFrames, Turing, CategoricalArrays, StatsBase, StatsPlots, Random


## read in data, it should have the following columns:
## flows, fromdist, todist, frompop, topop, distance, gdpcfrom, gdpcto, agegroup

Random.seed!(20240719)


ourdat = CSV.read("data/FlowData.csv",DataFrame)

ourdat.fromdist = categorical(ourdat.fromdist)
ourdat.todist = categorical(ourdat.todist)

ourdat.agegroup = categorical(ourdat.agegroup)
levels!(ourdat.agegroup,["unter18","18-25","25-30","30-50","50-65","über65"])


#meddist = median(ourdat.distance)

meddist = 293.0 ## median distance between districts according to Konstantin's printed report in km


### A Turing.jl sketch of a model

using Turing

@model function migration1(flows,fromdist,todist,frompop,topop,distance,gdpc,agegroup,Ndist,meddist)

    a ~ Gamma(5.0,1.0/4.0)
    b ~ Gamma(3.0,10.0/2.0)
    c ~ Gamma(5.0,1.0/4.0)
    d0 ~ Gamma(5.0,0.10/4.0)

    preds = frompop .* topop .* a ./ 1000.0 .* (1.0 .+ b ./ (dist ./ meddist .+ d0).^c)
    flows ~ arraydist([Poisson(p) for p in preds])
    
end


## sample the above model for 600 samples on 3 parallel chains with NUTS after doing 500 warmups and targeting 80% acceptance


fit1 = sample(migration1(ourdat.flows,levelcode.(ourdat.fromdist),levelcode.(ourdat.todist),ourdat.frompop,ourdat.topop,ourdat.distance,
    ourdat.gdpcfrom, ourdat.gdpcto, levelcode.(ourdat.agegroup),length(unique(ourdat.fromdist)),meddist),
    NUTS(500,.8), MCMCThreads(), 3, 600)

plot(fit1) |> display()

display(fit1)

### An alternative model is that we have a set of "desirability" numbers one set for each of the districts


function desir(d1,d2,d3,d4,d5,d6,age,from,to)
    if age == 1
        d1[to]/d1[from]
    elseif age == 2
        d2[to]/d2[from]
    elseif age == 3
        d3[to]/d3[from]
    elseif age == 4
        d4[to]/d4[from]
    elseif age == 5
        d5[to]/d5[from]
    elseif age == 6
        d6[to]/d6[from]                
    end
end

@model function migration2(flows,fromdist,todist,frompop,topop,distance,gdpcfrom,gdpcto,agegroup,Ndist,meddist,netactual)

    a ~ Gamma(5.0,1.0/4.0)
    b ~ Gamma(3.0,10.0/2.0)
    c ~ Gamma(5.0,1.0/4.0)
    d0 ~ Gamma(5.0,0.10/4.0)
    neterr ~ Exponential(0.1)

    desir1 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir2 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir3 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir4 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir5 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir6 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desires = [desir(desir1,desir2,desir3,desir4,desir5,desir6,
                        agegroup[i],fromdist[i],todist[i]) for i in 1:length(flows)]
    preds = frompop .* topop .* a ./ 1000.0 .* (1.0 .+ b ./ (dist ./ meddist .+ d0).^c) .* desires

    ins = zeros(typeof(preds[1]),Ndist)
    outs = zeros(typeof(preds[1]),Ndist)
    for n in 1:length(preds)
        ins[todist[i]] += preds[i]
        outs[fromdist[i]] -= preds[i]
    end
    netflows = ins-outs
    flows ~ arraydist([Poisson(p) for p in preds])
    netactual ~ arraydist(Normal.(netflows,neterr .* netflows))
end


fit2 = sample(migration2(ourdat.flows,levelcode.(ourdat.fromdist),levelcode.(ourdat.todist),ourdat.frompop,ourdat.topop,ourdat.distance,
    ourdat.gdpcfrom, ourdat.gdpcto, levelcode.(ourdat.agegroup), length(unique(ourdat.fromdist)),meddist),
    NUTS(500,.8), MCMCThreads(), 3, 600)

## too many parameters for this plot
#plot(fit2) |> display()

display(fit2)

mainparms2 = fit2[:,[:a,:b,:c,:d0],:]
plot(mainparms2) |> display()
