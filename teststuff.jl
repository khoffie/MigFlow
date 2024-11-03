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

function joint()
    lp = round(chain[:lp][end]; digits = 2)
    p = plot(chain[:X], chain[:Y], seriestype = :scatter,
             label = "LP = $(lp)")
    return p
end

function fit(;add = true)
    @model function foo()
        X ~ Normal(200, 10)
        Y ~ Normal(20, 10)
        if add 
            Turing.@addlogprob!(logpdf(Normal(1, .25), sqrt(X^2 + Y^2)))
        end
    end
    chain = Turing.sample(foo(), NUTS(), 1000)
    x = marginal(chain[:X])
    y = marginal(chain[:Y])
    j = joint()
    return j, x, y
end

function fitall()
    j, x, y = fit(; add = false)
    j2, x2, y2 = fit(; add = true)
    plt = plot(x, y, x2, y2, p2)
    display(plt)
    return p
end


j, x, y = fit(; add = false)
j2, x2, y2 = fit()

plot(j, j2)
plot(x, y, x2, y2, p2, layout = (3, 2))


