using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise

Revise.includet("./writeup/prior/jointprior.jl")

##        d = -.5(((x - μ₁) / σ₁)² + ((y - μ₂)/ σ₂)² + ((√(x² + y²) - μ₃) / σ₃)²)
        # z = ((√(x^2 + y^2) - μ₃) / σ₃)^2
        # d = -.5(((x - μ₁) / σ₁)^2 + ((y - μ₂)/ σ₂)^2 + ((z - μ₃) / σ₃ * √2)^2)

@model function custom_model()
    μ₁ = 10
    μ₂ = 10
    μ₃ = 1
    σ₁ = 10
    σ₂ = 10
    σ₃ = .01
    x ~ Normal(μ₁, σ₁)
    y ~ Normal(μ₂, σ₂)
##    z ~ Normal(μ₃, σ₃)
    function density(x, y)
        d = -.5(((x - μ₁) / σ₁)^2 + ((y - μ₂)/ σ₂)^2 + ((√(x^2 + y^2) - μ₃) / σ₃)^2)
        return d
    end
    Turing.@addlogprob!(density(x, y))
end
chain = Turing.sample(custom_model(), NUTS(), 1000)
plot(chain[:x], chain[:y], seriestype = :scatter)


