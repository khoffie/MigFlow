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

function joint(x, y)
    p = plot(x, y, seriestype = :scatter, label = "")
    return p
end

function allplots(chain)
    x = marginal(chain[:X])
    y = marginal(chain[:Y])
    z = marginal(chain[:Z])
    xy = joint(chain[:X], chain[:Y])
    xz = joint(chain[:X], chain[:Z])
    yz = joint(chain[:Y], chain[:Z])
    p = plot(x, y, z, xy, xz, yz, layout = (3, 3),
             title = ["X" "Y" "Z" "XY" "XZ" "YZ"])
    display(p)
    return p
end

function fit(; add = true)
    @model function foo()
        X ~ Normal(20, 10)
        Y ~ Normal(5, 5)
        Z = sqrt(X^2 + Y^2) 
        Z ~ Normal(10, 1)        
        if add 
            Turing.@addlogprob!(logpdf(Normal(1, 0.01), sqrt(X^2 + Y^2)))
        end
    end
    chain = Turing.sample(foo(), NUTS(), 1000)
    p = allplots(chain)
    return chain, p
end

chain1, p1 = fit(; add = false)
chain2, p2 = fit(; add = true)
savefig(p1, "./jointdist1.pdf")
savefig(p2, "./jointdist2.pdf")

