using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings

@model function custom_model()
    m1 = 2
    m2 = 20
    s1 = 1
    s2 = 1
    x ~ Normal(m1, s1)
    y ~ Normal(m2, s2)
    function custom_density(x, y)
        d = exp(-1/2((x / s1)^2 + (y / s2)^2 + ((sqrt(x^2 + y^2) -1) / .001)^2))
        return d
    end
    Turing.@addlogprob!(custom_density(x, y))
end
chain = Turing.sample(custom_model(), NUTS(), 1000)
plot(chain[:x], chain[:y], seriestype = :scatter)


@model function custom_model()
    function custom_density(x, y)
##        d = exp(-1/2((x/10)^2))
        d = exp(-.5 * (x^2 + y^2 - 1) / .01)
        return d
    end
    x ~ Normal(20, 1)
    y ~ Normal(20, 10)
    Turing.@addlogprob!(custom_density(x, y))
end

marginal(chain[:x])
marginal(chain[:y])

