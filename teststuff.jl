using Turing: @addlogprob!
using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings

function marginal(x)
    m = round(mean(x); digits = 2)
    sd = round(std(x); digits = 2)
    p = density(x, label = L"\mu = %$m, \sigma = %$sd")
    display(p)
    return p
end

function joint(chain)
    lp = round(chain[:lp][end]; digits = 2)
    p = plot(chain[:X], chain[:Y], seriestype = :scatter,
             label = "LP = $(lp)")
    return p
end

function fit(; add = true)
    @model function foo()
        X ~ Normal(20, 10)
        Y ~ Normal(5, 5)
        if add 
            Turing.@addlogprob!(logpdf(Normal(1, 0.01), sqrt(X^2 + Y^2)))
        end
    end
    chain = Turing.sample(foo(), NUTS(), 1000)
    return chain
end

function fitall()
    chain1 = fit(; add = false)
    chain2 = fit(; add = true)
    x1 = marginal(chain1[:X])
    y1 = marginal(chain1[:Y])
    x2 = marginal(chain2[:X])
    y2 = marginal(chain2[:Y])
    j = joint(chain2)
    p = plot(x1, y1, x2, y2, j, layout = (3, 2), title = ["X" "Y" "X" "Y" "XY"])
    return p
end

p = fitall()

