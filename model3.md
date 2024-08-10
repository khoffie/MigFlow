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

