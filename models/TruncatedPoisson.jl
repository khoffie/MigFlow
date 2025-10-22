struct TruncatedPoisson{T<:Real} <: DiscreteUnivariateDistribution λ::T end

function Distributions.logpdf(d::TruncatedPoisson, k::Real)
    if k < 1
        return -Inf
    else
        return logpdf(Poisson(d.λ), k) - log1mexp(-d.λ)
    end
end

function Distributions.rand(rng::AbstractRNG, d::TruncatedPoisson)
    k = 0
    while k < 1
        k = rand(rng, Poisson(d.λ))
    end
    return k
end

# not mathematical support 1 to Inf but reasonable range for plotting
Distributions.support(d::TruncatedPoisson) = 1:ceil(Int, d.λ + 5√(d.λ))
Distributions.minimum(::TruncatedPoisson) = 1
Distributions.maximum(::TruncatedPoisson) = Inf

function Statistics.quantile(d::TruncatedPoisson, p::Real)
    # used by Makie
    λ = d.λ
    P0 = pdf(Poisson(λ), 0)
    # Transform back to untruncated Poisson quantile
    q_pois = quantile(Poisson(λ), p * (1 - P0) + P0)
    return max(q_pois, 1)
end
