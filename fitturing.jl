
using CSV,DataFrames,Turing, CategoricalArrays


## read in data, it should have the following columns:
## flowcount, fromdist, todist, frompop, topop, distance, gdppcfrom,gdppcto, agegroup


ourdat = CSV.read("data/FlowData.csv",DataFrame)

ourdat.

### A Turing.jl sketch of a model

using Turing

@model migration1(flows,fromdist,todist,frompop,topop,distance,gdpc,agegroup,Ndist,meddist)

    a ~ Gamma(5.0,1.0/4.0)
    b ~ Gamma(5.0,1.0/4.0)
    c ~ Gamma(5.0,1.0/4.0)
    d0 ~ Gamma(5.0,1.0/4.0)

    preds = [frompop * topop * a/1000.0 * (1.0 + b/(dist/meddist + d0)^c) for _ in flows]
    flows ~ arraydist([Poisson(p) for p in preds])
    
end

function desir(d1,d2,d3,d4,age,from,to)
    if age == 1
        d1[to]/d1[from]
    elseif age == 2
        d2[to]/d2[from]
    elseif age == 3
        d3[to]/d3[from]
    elseif age == 4
        d4[to]/d4[from]
    end
end

@model migration2(flows,fromdist,todist,frompop,topop,distance,gdpcfrom,gdpcto,agegroup,Ndist,meddist)

    a ~ Gamma(5.0,1.0/4.0)
    b ~ Gamma(5.0,1.0/4.0)
    c ~ Gamma(5.0,1.0/4.0)
    d0 ~ Gamma(5.0,1.0/4.0)
    desir1 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir2 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir3 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir4 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desires = [desir(desir1,desir2,desir3,desir4,agegroup[i],fromdist[i],todist[i]) for i in 1:length(flows)]
    preds = [frompop * topop * a/1000.0 * (1.0 + b/(dist/meddist + d0)^c) * desires ] 
    flows ~ arraydist([Poisson(p) for p in preds])
end

