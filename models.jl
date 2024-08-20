
#=

migration1 uses a symmetric formulation and has only distance function

=#


@model function migration1(flows, fromdist, todist, frompop, topop, distance,
    agegroup, Nages, Ndist, meddist)

    a ~ filldist(Gamma(5.0, 10.0/4.0), Nages)
    b ~ filldist(Gamma(3.0, 1.0/2.0), Nages)
    c ~ filldist(Gamma(5.0, 1.0/4.0), Nages)
    d0 ~ filldist(Gamma(5.0, 0.10/4.0), Nages)

    preds = [frompop[i] * topop[i] * a[agegroup[i]] / 1000.0 * (1.0 + b[agegroup[i]] / (distance[i] / meddist + d0[agegroup[i]])^c[agegroup[i]]) for i in 1:length(flows)]

    if typeof(a[1]) != Float64
        @printf "a1 = %.2f, b1 = %.2f, c1 = %.2f, d0_1 = %.2f\n"  a[1].value b[1].value c[1].value d0[1].value
    end
    if any(isnan,preds)
        println("NaN in predictions")
    end
        flows ~ arraydist([Poisson(p) for p in preds])
    return preds
end





#=

migration2 needs inputs:


=#




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






#=


Model 3: migration3 needs as inputs:

etc etc


To see the distance function: 

https://www.desmos.com/calculator/jhrgbmw9dd
=#


@model function migration3(flows, allmoves, fromdist, todist,
    frompop, topop, popgerm, distance,
    agegroup, Nages,
    xcoord, ycoord,density,distpop,
    Ndist, meddist, netactual, ncoefs)

    a ~ filldist(Normal(0.0,2.75),Nages) # center a around a typical value found by plotting the data and adjusting until the overall shift on log scale is about right
    c ~ filldist(Gamma(10.0, 1.5/9.0),Nages) # c decay rate parameter between about 1 and 4 ish
    d0 ~ filldist(Gamma(5.0, 2.0/4.0),Nages)
    dscale ~ filldist(Gamma(8.0,0.5/7.0),Nages) ## sensitive falls to same size as insensitive somewhere around 0.5 times the median distance between districts
    neterr ~ Gamma(3.0, 5/2.0) ## this is in percent
    kd ~ MvNormal(fill(0.0, Nages), (log(15.0) / 0.5) / 2 * ones(Nages)) # density ranges mostly in the range -0.5 to 0.5, so a full-scale change in density could multiply the flow by around 5.0

    ## priors for chebychev polys parameters
    desirecoefs ~ MvNormal(zeros(ncoefs*Nages), 2.0 .* ones(ncoefs*Nages)) ## new scaling we are dividing by 10 below
    desirecoefsre = reshape(desirecoefs, (ncoefs,Nages)) ## vec -> matrix
    ### creates functions from x and y coords, Fun callable chebychev polynomial, kind of a function
    ### array of funs per age group
    desfuns = [Fun(Chebyshev(300.0 .. 1000.0) * Chebyshev(5000.0 .. 6200.0), desirecoefsre[:,i] ./ 10) for i in 1:Nages]
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

    ## indiviudal flows, the 4.5 is a number we got by approximately centering the a values to make them more interpretable
    distscale = dscale .* meddist

    preds = [frompop[i] * logistic( -3.0 + log(topop[i] / popgerm) + a[agegroup[i]] +
                log1p(1.0 / (distance[i] / distscale[agegroup[i]] + d0[agegroup[i]]/100.0)^c[agegroup[i]]) + desires[i])
                    for i in 1:length(flows)]

    if typeof(a[1]) != Float64
        @printf "a1 = %.2f, c1 = %.2f, d0_1 = %.2f, neterr = %.2f\n"  a[1].value c[1].value d0[1].value neterr.value
    end


    #if typeof(a[1]) != Float64
    #    @printf "a1 = %.2f, b1 = %.2f, c1 = %.2f, d0_1 = %.2f, neterr = %.2f\n"  a[1] b[1] c[1] d0[1] neterr
    #end 

    if any(isnan,preds)
        println("NaN in predictions")
    end
    ## matrix dist, age
    netflows = calcnet(preds,fromdist,todist,agegroup,Nages,Ndist)

    predmoves = sum(preds) # total predicted flow
    allmoves ~ Normal(predmoves, .01 * predmoves) ## our total move predictions should be about right by around 1%

    ## flows ~ poisson(expectation)
    flows ~ arraydist([Poisson(p) for p in preds])

    ### both netactual and flows we optimize together. neterr   
    ### determinies the importance of netactual for optimization. the
    ### "best prediction" for flows is biased by netactual.
    ### bathtube

    ## bubble under rug, how do I get the bubble to a certain spot by
    ## stomping on the rug?
    neterrfrac = neterr/1000 ## we rescaled it to tenths of a percent of district population
    netactual ~ arraydist([Normal(netflows[dist,age],neterrfrac*distpop[dist]) for dist in 1:Ndist, age in 1:Nages]) 
    return((preds,netflows))
end
