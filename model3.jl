@model function migration3(flows, allmoves, fromdist, todist,
                           frompop, topop, popgerm, distance,
                           agegroup, Nages,
                           xcoord, ycoord,density,distpop,
                           Ndist, meddist, netactual, ncoefs)
    a ~ filldist(Normal(0.0,1.0),Nages)
    b ~ filldist(Gamma(3.0, 1.0/2.0),Nages)
    c ~ filldist(Gamma(5.0, 2.0/4.0),Nages)
    d0 ~ filldist(Gamma(5.0, 0.2/4.0),Nages)
    neterr ~ Gamma(3.0, .05/2.0)
    logisticconst ~ Uniform(-12.0,12.0) # This constant isn't easy to figure out because log(topop[i]/popgerm) is numbers in the range maybe -10 to -4 
    kd ~ MvNormal(fill(0.0, Nages), (log(5.0) / 0.5) / 2 * ones(Nages)) # density ranges mostly in the range -0.5 to 0.5, so a full-scale change in density could multiply the flow by around 5.0

    ## priors for chebychev polys parameters
    desirecoefs ~ MvNormal(zeros(ncoefs*Nages), 1.0 .* ones(ncoefs*Nages))
    desirecoefsre = reshape(desirecoefs, (ncoefs,Nages)) ## vec -> matrix
    ### creates functions from x and y coords, Fun callable chebychev polynomial, kind of a function
    ### array of funs per age group
    desfuns = [Fun(Chebyshev(300.0 .. 1000.0) * Chebyshev(5000.0 .. 6200.0), desirecoefsre[:,i]) for i in 1:Nages]
    desvals = [desfuns[age](xcoord,ycoord) for (xcoord,ycoord) in zip(xcoord,ycoord), age in 1:Nages]

    ## as prior for the coefficients, let's tell it that desvals should be near 0 and have standard deviation 1 ish.
    Turing.@addlogprob!(logpdf(Normal(0.0,0.25),mean(desvals)))
    Turing.@addlogprob!(logpdf(Exponential(1.0),std(desvals)))

    ## desirability for given age at coordinates as ratio of dest / from
    desires = [(kd[agegroup[i]] * density[todist[i]] +
                    desvals[todist[i], agegroup[i]]) -
                    (kd[agegroup[i]] * density[fromdist[i]] +
                        desvals[fromdist[i], agegroup[i]])
               for i in 1:length(flows)]
    ## indiviudal flows
    preds = [frompop[i] * logistic(logisticconst + log(topop[i] / popgerm) + a[agegroup[i]] +
        log1p(b[agegroup[i]] / (distance[i] / meddist + d0[agegroup[i]])^c[agegroup[i]]) + desires[i])
             for i in 1:length(flows)]

    if typeof(a[1]) == Float64
        @printf "a1 = %.2f, b1 = %.2f, c1 = %.2f, d0_1 = %.2f, neterr = %.2f\n"  a[1] b[1] c[1] d0[1] neterr
    end
    if typeof(a[1]) != Float64
        @printf "a1 = %.2f, b1 = %.2f, c1 = %.2f, d0_1 = %.2f, neterr = %.2f\n"  a[1].value b[1].value c[1].value d0[1].value neterr.value
    end
    
    if any(isnan,preds)
        println("NaN in predictions")
    end
    ## matrix dist, age
    netflows = calcnet(preds,fromdist,todist,agegroup,Nages,Ndist)

    predmoves = sum(preds) # total predicted flow
    allmoves ~ Normal(predmoves, .01*predmoves) ## our total move predictions should be about right by around 1%

    # total net flow as fraction of germany should be very close to zero: this is part of the overall prior on all parameters
    # we might want to make the scale here be a parameter but let's start with an allowable imbalance around 50 people / million people in Germany
    Turing.@addlogprob!(logpdf(Normal(0.0,50.0/1e6), sum(netflows) / popgerm))

    ## flows ~ poisson(expectation)
    flows ~ arraydist([Poisson(p) for p in preds])
    
    ### both netactual and flows we optimize together. neterr
    ### determinies the importance of netactual for optimization. the
    ### "best prediction" for flows is biased by netactual.
    ### bathtube

    ## bubble under rug, how do I get the bubble to a certain spot by
    ## stomping on the rug?
    netactual ~ arraydist([Normal(netflows[dist,age],neterr*distpop[dist]) for dist in 1:Ndist, age in 1:Nages])
    return((preds,netflows))
end
