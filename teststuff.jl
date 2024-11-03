using Turing: @addlogprob!
lusing Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots


function fit(;add = true)
    @model function foo()
        X ~ Normal(200, 1)
        Y ~ Normal(200, 1)
        if add 
            Turing.@addlogprob!(logpdf(Normal(1, .25), sqrt(X^2 + Y^2)))
        end
    end
    chain = Turing.sample(foo(), NUTS(), 1000)
    lp = round(chain[:lp][end]; digits = 2)
    p = plot(chain[:X], chain[:Y], seriestype = :scatter,
             label = "LP = $(lp)")
    display(p)
    return chain
end

chain = fit(add = true)




