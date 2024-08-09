@model function migration3(flows, fromdist, todist,
                           frompop, topop, popgerm, distance,
                           agegroup, Nages,
                           xcoord, ycoord,density,
                           Ndist, meddist, netactual, ncoefs)
    a ~ filldist(Normal(0.0,4.0),Nages)
    b ~ filldist(Gamma(3.0, 1.0/2.0),Nages)
    c ~ filldist(Gamma(5.0, 1.0/4.0),Nages)
    d0 ~ filldist(Gamma(5.0, 1.0/4.0),Nages)
    neterr ~ Gamma(3.0, 2.0/2.0)
    kd ~ MvNormal(fill(0.0,Nages),0.1*ones(Nages))
    logisticconst ~ Normal(-4.0,2.0) # logistic(-4.0) ~ 0.017 flows are typically on order 5% or less

    ## priors for chebychev polys parameters
    desirecoefs ~ MvNormal(zeros(ncoefs*Nages), 1.0 .* ones(ncoefs*Nages))  
    desirecoefsre = reshape(desirecoefs, (ncoefs,Nages)) ## vec -> matrix
    ### creates functions from x and y coords, Fun callable chebychev polynomial, kind of a function
    ### array of funs per age group
    desfuns = [Fun(Chebyshev(200000.0 .. 600000.0) * Chebyshev(5e6 .. 6.2e6), desirecoefsre[:,i]) for i in 1:Nages]
    desvals = [desfuns[age](xcoord,ycoord) for (xcoord,ycoord) in zip(xcoord,ycoord), age in 1:Nages]
    ## desirability for given age at coordinates as ratio of dest / from
    desires = [(kd[agegroup[i]]*density[todist[i]] +
                    desvals[todist[i],agegroup[i]]) - 
                    (kd[agegroup[i]]*density[fromdist[i]] +
                        desvals[fromdist[i],agegroup[i]])
               for i in 1:length(flows)]
    ## indiviudal flows
    preds = [frompop[i] * logistic(logisticconst + log(topop[i] / popgerm) + a[agegroup[i]] +
        log1p(b[agegroup[i]] / (distance[i] / meddist + d0[agegroup[i]]) * c[agegroup[i]]) + desires[i])
             for i in 1:length(flows)]

    if typeof(a[1]) != Float64
        @printf "a1 = %.2f, b1 = %.2f, c1 = %.2f, d0_1 = %.2f, neterr = %.2f\n"  a[1].value b[1].value c[1].value d0[1].value neterr.value
    end
    
    if any(isnan,preds)
        println("NaN in predictions")
    end
    ## matrix dist, age
    netflows = calcnet(preds,fromdist,todist,agegroup,Nages,Ndist)

    # total net flow as fraction of germany should be very close to zero: this is part of the overall prior on all parameters
    Turing.@addlogprob!(logpdf(Normal(0.0,0.005),sum(netflows)/popgerm))

    ## flows ~ poisson(expectation)
    flows ~ arraydist([Poisson(p) for p in preds])
    # for c in axes(netactual,2)
    #     netactual[:,c] ~ MvNormal(netflows[:,c],neterr .* abs.(netflows[:,c]))
    # end
    
    ### both netactual and flows we optimize together. neterr
    ### determinies the importance of netactual for optimization. the
    ### "best prediction" for flows is biased by netactual.
    ### bathtube

    ## bubble under rug, how do I get the bubble to a certain spot by
    ## stomping on the rug?
    netactual ~ arraydist(Normal.(netflows, neterr .* abs.(netflows)))
    return((preds,netflows))
end
