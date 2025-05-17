function distnorm(data::NamedTuple)
    df        = data.df
    districts = data.districts
    Ndist = length(districts)
    dist = data.distances
    dffull    = data.dffull
    Y       = df.flows
    from    = levelcode.(categorical(df.fromdist))
    to      = levelcode.(categorical(df.todist))
    AP      = unique(df, :fromdist).frompop ## population of agegroup
    poporig = districts.pop
    P       = log.(districts.pop ./ 153000) # median pop
    D       = df.dist ./ 100.0
    N       = length(Y)
    Nfull   = length(dffull.fromdist)
    fromfull = levelcode.(categorical(dffull.fromdist))
    tofull   = levelcode.(categorical(dffull.todist))
    Dfull    = dffull.dist
    nfrom    = length(unique(from))
    radius   = sqrt.(districts.pop ./ districts.density ./ 2π)
    data = (; Y, D, from, to, AP, P, poporig)

    @model function model(Y, from, to, D, AP, P, N, Nfull, fromfull, tofull,
                          Dfull, nfrom, radius)
        α ~ Normal(-8, 1)
        β ~ Gamma(1, 1)
        γ_raw ~ Gamma(15, .2); γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1); ϕ = ϕ_raw / 100
        δ_raw ~ Gamma(10, 1); δ = δ_raw / 100

        T = eltype(γ)  # to get dual data type for AD
        denom = zeros(T, nfrom)
        ps = Vector{T}(undef, N)

        desmat = [exp(desirability(P[to[i]], dist[i,j], ϕ, δ, γ))
                  for i in 1:Ndist, j in 1:Ndist]
        denom = [sum(desmat[i,j] for i in 1:Ndist) for j in 1:Ndist]

        @inbounds for i in 1:N
            ps[i] = AP[from[i]] * exp(α) * desmat[to[i], from[i]] / denom[from[i]]
        end
        Y .~ Poisson.(ps)
        return ps
    end

    mdl = model(Y, from, to, D, AP, P, N, Nfull, fromfull,
                tofull, Dfull, nfrom, radius)
    lb = [-20, -100, 10, 0, 1]
    ub = [0, 10, 100, 99, 100]

    return (; mdl, lb, ub, data)
end

desirability(P, D, ϕ, δ, γ) = P + log((ϕ + (1 - ϕ) / (D + δ)^γ))
