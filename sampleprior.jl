@model function usprior(ndenscoef, ncoefs, xcoord, ycoord)
    xmin = minimum(xcoord)
    xmax = maximum(xcoord)
    ymin = minimum(ycoord)
    ymax = maximum(ycoord)
    a ~ Normal(-14.0, 7)
    c ~ Gamma(10.0, 1.5 / 9.0)
    d0 ~ Gamma(5.0, 2.0 / 4.0)
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
    return desvals
end

districts = CSV.read("./data/districts.csv", DataFrame)
uschain = Turing.sample(usprior(36, 36, districts.xcoord, districts.ycoord), NUTS(), 10)

function sampleprior(model, values, params)
    n = 10
    df = DataFrame(Matrix{Float64}(undef, 401, 10), :auto)
    for i in 1 : n
        df[:, i] = generated_quantities(model, values[i, 1 : 77, 1], params)
    end
    return df
end

samples = sampleprior(
    usprior(36, 36, districts.xcoord, districts.ycoord),
    uschain.value.data,
    uschain.name_map.parameters)

samples[:, "distcode"] = districts[:, "distcode"]
CSV.write("./manuscript_input/priorsamples.csv", samples)
