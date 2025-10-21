struct TruncatedPoisson{T<:Real} <: DiscreteUnivariateDistribution λ::T end

function Distributions.logpdf(d::TruncatedPoisson, k::Real)
    if k < 1
        return -Inf
    else
        λ = d.λ
        ## return k * log(λ) - λ - logfactorial(k) - log1mexp(-λ)
        return logpdf(Poisson(λ), k) - log1mexp(-λ)
    end
end

function Distributions.rand(rng::AbstractRNG, d::TruncatedPoisson)
    k = 0
    while k < 1
        k = rand(rng, Poisson(d.λ))
    end
    return k
end
