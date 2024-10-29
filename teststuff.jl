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
