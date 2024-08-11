# Migration3 model

The migration3 Turing model works as follows (as of 2024/9/9)

For each flow from location A to location B we have a prediction. This prediction is in terms of the fraction of location A which leaves for location B

So the number of people going A -> B would be the population of A times this fraction:

```{julia}

preds = [frompop[i] * logistic(logisticconst + log(topop[i] / popgerm) + a[agegroup[i]] +
        log1p(b[agegroup[i]] / (distance[i] / meddist + d0[agegroup[i]]) * c[agegroup[i]]) + desires[i])
             for i in 1:length(flows)]

```

the logistic function is $f(x) = 1/(1+\exp(-x))$ except implemented in the computer in a numerically stable way.

The logistic function behaves like exp(x) when x is less than about -4, but saturates at 1 for large x

We expect almost all districts in Germany have a total outflow to any one location less than a couple percent of their current population. so we think in the logistic curve logistic(x) should have x in the range less than -3 (where logistic(-3) ~ 0.05)

Given that, we start with a constant value "logisticconst" which has a prior Normal(-4.0,2.0). 

Then each age group has an a value which adds to this, allowing variation in overall migrativity magnitude per age group, and then there's a function of distance... and then there's the difference between chebyshev polynomials etc expressed in "desires[i]"

In addition to all these requirements, we require that the total migration in germany should be about 0 (because people leaving one place are going to another place). We express this as the total netflow as a fraction of the total german population is very close to 0.

```{julia}
    Turing.@addlogprob!(logpdf(Normal(0.0,0.005),sum(netflows)/popgerm))
```

This term in the log probability density can be thought of as a dependent part of a prior, depending on all the parameters together, such as the chebyshev, a,b,c, etc parameters. They should all combine in such a way that the total migration is 0, no overall imbalances. We don't require this to be precisely 0 but strongly downweight any result that has more than about 0.5% error

This is a critical restriction on the form of the function!

Our likelihood is based on two kinds of observations:

1) observed individual flows
```{julia}
    flows ~ arraydist([Poisson(p) for p in preds])
```

2) observed net flows
```{julia}
    netactual ~ arraydist(Normal.(netflows, neterr .* abs.(netflows)))
```

## To Do and Notes

for 2024-08-12

We have converted neterr to a percentage scale, and now we need the district population
```{julia}
    neterrfrac = neterr/100 ## we rescaled it to percent
    netactual ~ arraydist([Normal(netflows[dist,age],neterrfrac*distpop[dist]) for dist in 1:Ndist, age in 1:Nages]) # neterr is in percent now
```

I've also constrained the Chebyshev coefficients more, first I relaxed the prior a little:

```{julia}
    desirecoefs ~ MvNormal(zeros(ncoefs*Nages), 2.0 .* ones(ncoefs*Nages)) ## new scaling we are dividing by 10 below

```
But then I rescaled the coefficients by a factor of 10 before using them:

```{julia}
    desfuns = [Fun(Chebyshev(300.0 .. 1000.0) * Chebyshev(5000.0 .. 6200.0), desirecoefsre[:,i] ./ 10) for i in 1:Nages]

```

So now we'll want to understand that these coefficients are 1/10 what they used to be, so now we need to consider that in the lower and upper bounds, they should be around 10x what they were.

We also needed to add the district population to the model, which I tried to do as `dists.pop`

```
    model3 = migration3(dt2.flows, sum(dt2.flows), levelcode.(dt2.fromdist), levelcode.(dt2.todist),
                        dt2.frompop, dt2.topop, popgerm, dt2.distance,
                        levelcode.(dt2.agegroup),
                        Nages,
                        dists.xcoord, dists.ycoord, distdens,dists.pop,
                        Ndist, meddist, netactual, ncoefs)

```

We may want to expand the neterr parameter into one per age group, but for now let's see if this helps.

So, given this new set of fixes I think we should try as follows:

1) Run some short optimizations using kd and cheby coefficients near 0. find an optimum and make graphs of log(actual/pred) vs distance and log(actual) vs log(pred) for pred and actual > 3

2) Relax the constraints on kd to something like +- 1.5 and keep cheby +- 0.2 and run optimization and then do the same graphs as above. Plot an additional single map with residual values.

3) Relax chebyshev to +- 0.5 and 1.0 and 2.0 at each one time plot maps and fits. Run longer for these guys.

4) Find a set of initial conditions that make sense from (3) and use them as inits for a sample, NUTS(500,0.8) and 100 real samples, using MCMCThreads() and 3 threads. This may take a long time

CONTINGENCY: If things are problematic, we may need to tune the error size, and/or 