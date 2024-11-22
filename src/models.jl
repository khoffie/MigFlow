@model function usmodel(flows, allmoves, fromdist, todist, medcpop, distance,
    xcoord, xmin, xmax, ycoord, ymin, ymax, logreldens,densmin, densmax, distpop,
    Ndist, meddist, ncoefs, ndenscoef, positive_only)

    a ~ Normal(-14.0, 7)
    c ~ Gamma(10.0, 1.5 / 9.0)
    d0 ~ Gamma(5.0, 2.0 / 4.0)
    e ~ Gamma(5.0, 2.0 / 4.0)
    dscale ~ Gamma(20.0, 5.0 / 19.0)
    ktopop ~ Normal(0.0, 5.0) # remove? model should be linear in topop
    kd ~ MvNormal(zeros(ndenscoef), fill(40.0, ndenscoef))
    desirecoefs ~ MvNormal(zeros(ncoefs), fill(50.0, ncoefs))

    ## The operator * creates a 2D space from both 1D spaces as
    ## defined by Chebyshev(xmin .. xmax) and Chebyshev(ymin
    ## .. ymax). At each point (x, y) in the 2D plane, we have a
    ## product between two Cheby polynomials. 
    desfun = Fun(ApproxFun.Chebyshev(xmin .. xmax) * ApproxFun.Chebyshev(ymin .. ymax), desirecoefs ./ 10)
    ## evalueates the 2D desirability cheby for each district
    desvals = [desfun(xcoord, ycoord) for (xcoord, ycoord) in zip(xcoord, ycoord)]

    mu = fill(0.0, length(desvals))
    sigma = fill(3.0, length(desvals))
    ## logpdf returns the log of density of MvN at all 401 districts
    Turing.@addlogprob!(logpdf(MvNormal(mu, sigma), desvals))
    ## tamp down the corners of the map
    Turing.@addlogprob!(logpdf(MvNormal(fill(0.0, 4), fill(0.125, 4)),
                               [desfun(xmin, ymin),
                                desfun(xmin, ymax),
                                desfun(xmax, ymin),
                                desfun(xmax, ymax)]))

    densfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax), kd ./ 10)
    distscale = dscale .* meddist
    preds = [distpop[fromdist[i]] *
        logistic(densfun(logreldens[fromdist[i]], logreldens[todist[i]]) +
        a + (ktopop / 10.0) *
        log(distpop[todist[i]] / medcpop + e / 100.0) +
        log1p(1.0 / (distance[i] / distscale + d0 / 100.0) ^ c) +
        desvals[todist[i]] - desvals[fromdist[i]]) for i in eachindex(flows)]
    #logpreds = log.(preds)
    ## the logarithm of predicted rates is in the range of about -9 to +10 this is more or less a prior on 
    ## the kd and desfun parameters
#    Turing.@addlogprob!(logpdf(MvNormal(fill(1.0,length(logpreds)),fill(2.0,length(logpreds))),logpreds)) 

    sumpreds = sum(preds)
    allmoves ~ Normal(sumpreds, 0.01 * sumpreds)
    if any(isnan, preds)
        println("NaN in predictions")
    end
    if positive_only
        flows ~ arraydist([truncated(Poisson(p), 1, Inf) for p in preds])        
    else
        flows ~ arraydist([Poisson(p) for p in preds])
    end
    return (preds = preds)
end

function germmodel(flows, districts, positive_only)
    model = usmodel(flows.flows, sum(flows.flows),
                   levelcode.(flows.fromdist), levelcode.(flows.todist),
                   median(districts.pop), flows.dist,
                   districts.xcoord, minimum(districts.xcoord), maximum(districts.xcoord),
                   districts.ycoord, minimum(districts.ycoord), maximum(districts.ycoord),
                   districts.logreldens, minimum(districts.logreldens), maximum(districts.logreldens),
                   districts.pop, nrow(districts),
                    100.0, 36, 36, pos_only) ## nothing == netactual, we're not using it anymorec
    return model
end
