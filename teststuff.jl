using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise

Revise.includet("./writeup/prior/jointprior.jl")

##        d = -.5(((x - μ₁) / σ₁)² + ((y - μ₂)/ σ₂)² + ((√(x² + y²) - μ₃) / σ₃)²)
        # z = ((√(x^2 + y^2) - μ₃) / σ₃)^2
        # d = -.5(((x - μ₁) / σ₁)^2 + ((y - μ₂)/ σ₂)^2 + ((z - μ₃) / σ₃ * √2)^2)

function fitmodel(add)
    chain = Turing.sample(custom_model(add), NUTS(), 1000)
    xy = plot(chain[:x], chain[:y], seriestype = :scatter, label = L"(X, Y)")
    x = density(chain[:x], label = L"X")
    y = density(chain[:y], label = L"X")
    xyd = density(sqrt.(chain[:x].^2 .+ chain[:y].^2), label = L"\sqrt{(X^2 + Y^2)}")
    p = plot(x, y, xy, xyd)
    display(p)
    return chain, p
end

@model function custom_model(add = true)
    μ₁ = 0
    μ₂ = 0
    μ₃ = 1
    σ₁ = 10
    σ₂ = 10
    σ₃ = 0.01
    x ~ Normal(μ₁, σ₁)
    y ~ Normal(μ₂, σ₂)
##    z ~ Normal(μ₃, σ₃)
    function density(x, y)
        d = -.5(((x - μ₁) / σ₁)^2 + ((y - μ₂)/ σ₂)^2 + ((√(x^2 + y^2) - μ₃) / σ₃)^2)
        return d
    end
    if add; Turing.@addlogprob!(density(x, y)); end
end

fitmodel(true)
