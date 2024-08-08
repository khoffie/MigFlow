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

@model function migration2(flows,fromdist,todist,frompop,topop,distance,agegroup,Nages,Ndist,meddist, netactual)
    a ~ filldist(Gamma(5.0, 10.0/4.0),Nages)
    b ~ filldist(Gamma(3.0, 1.0/2.0),Nages)
    c ~ filldist(Gamma(5.0, 1.0/4.0),Nages)
    d0 ~ filldist(Gamma(5.0, 0.10/4.0),Nages)
    neterr ~ Gamma(3.0, 2.0/2.0)

    desir ~ filldist(Gamma(400.0, 100.0/399.0), Ndist, Nages)

    desires = [desir[todist[i],agegroup[i]]/desir[fromdist[i],agegroup[i]] for i in 1:length(flows)]
    preds = [frompop[i] * topop[i] * a[agegroup[i]] / 1000.0 * (1.0 + b[agegroup[i]] / (distance[i] / meddist + d0[agegroup[i]])^c[agegroup[i]]) * desires[i] for i in 1:length(flows)]

    if typeof(a[1]) != Float64
        @show a[1].value,b[1].value,c[1].value,d0[1].value,neterr.value
    end
    if any(isnan,preds)
        println("NaN in predictions")
    end
    netflows = calcnet(preds,fromdist,todist,agegroup,Nages,Ndist)
    flows ~ arraydist([Poisson(p) for p in preds])
    netactual ~ arraydist(Normal.(netflows,neterr .* abs.(netflows)))
##    netactual ~ MvNormal(netflows,neterr .* abs.(netflows))
    return((preds,netflows))
end
