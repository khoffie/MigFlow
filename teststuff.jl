using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots
germ = loadallGermData(; sample = false)
us = loadallUSdata()


parnames = [["a", "c", "d0", "dscale", "ktopop"];
            ["kd[$i]" for i in 1:36];
            ["desirecoefs[$i]" for i in 1:36]]

lb = [[-60.0, 0.0, 0.0, 1.0, -10.0];
      fill(-50.50, 36);
      fill(-50.50, 36)]

ub = [[60.0, 20.0, 10.0, 15.0, 10.0];
      fill(50.50, 36);
      fill(50.50, 36)]

ini = rand(Normal(0.0, 0.10), length(ub))
ini[1:5] .= [-7.6, 1.81, 1.5, 5.0, 3.5]

new = DataFrame([ :pars => parnames, :lb => lb, :inits => ini, :ub => ub])


repeat(inits, 1, 4)
test = Iterators.repeated(inits.inits, 4)


 mhsamp = Turing.sample(alldata.model, NUTS(), MCMCThreads(),
                               nsamples, nchains, thinning = thinning,
                               initial_params = Iterators.repeated(vals.optis),
                               verbose = true, progress = true)

# Define a simple Normal model with unknown mean and variance.
@model function gdemo(x, y)
    s² ~ InverseGamma(2, 3)
    m ~ Normal(0, sqrt(s²))
    x ~ Normal(m, sqrt(s²))
    return y ~ Normal(m, sqrt(s²))
end


inits = [1, 2, 1, 2]
p1 = Turing.sample(gdemo(missing, missing), NUTS(), 100, thinning = 5,
                   initial_params = inits)
p1 = Turing.sample(gdemo(missing, missing), NUTS(),  MCMCThreads(), 100, 2, thinning = 5,
                   initial_params = Iterators.repeated(inits))

