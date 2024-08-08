@model function migration3(flows, fromdist, todist,
                           frompop, topop, distance,
                           agegroup, Nages,
                           xcoord, ycoord,
                           Ndist, meddist, netactual, ncoefs)
    a ~ filldist(Gamma(5.0, 10.0/4.0),Nages)
    b ~ filldist(Gamma(3.0, 1.0/2.0),Nages)
    c ~ filldist(Gamma(5.0, 1.0/4.0),Nages)
    d0 ~ filldist(Gamma(5.0, 0.10/4.0),Nages)
    neterr ~ Gamma(3.0, 2.0/2.0)

    
    desirecoefs ~ MvNormal(zeros(ncoefs),ones(ncoefs)) 
    desfun = Fun(Chebyshev(200000.0 .. 600000.0) * Chebyshev(5e6 .. 6.2e6), desirecoefs)
    
    desires = [exp(desfun(xcoord[todist[i]],ycoord[todist[i]]) -
        desfun(xcoord[fromdist[i]],ycoord[fromdist[i]]))
               for i in 1:length(flows)]
    preds = [frompop[i] * topop[i] * a[agegroup[i]] / 1000.0 *
        (1.0 + b[agegroup[i]] / (distance[i] / meddist + d0[agegroup[i]])^c[agegroup[i]]) * desires[i]
             for i in 1:length(flows)]

    if typeof(a[1]) != Float64
        @printf "a1 = %.2f, b1 = %.2f, c1 = %.2f, d0_1 = %.2f, neterr = %.2f"  a[1].value b[1].value c[1].value d0[1].value neterr.value
    end
    
    if any(isnan,preds)
        println("NaN in predictions")
    end
    netflows = calcnet(preds,fromdist,todist,agegroup,Nages,Ndist)
    flows ~ arraydist([Poisson(p) for p in preds])
    # for c in axes(netactual,2)
    #     netactual[:,c] ~ MvNormal(netflows[:,c],neterr .* abs.(netflows[:,c]))
    # end
     netactual ~ arraydist(Normal.(netflows,neterr .* abs.(netflows)))
    return((preds,netflows))
end
