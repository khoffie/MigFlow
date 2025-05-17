 function dan(data::NamedTuple; scaleo)
    df        = data.df
    districts = data.districts
    dffull    = data.dffull
    Ndist   = length(districts.distcode)
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
    DM      = make_distmat(districts.distcode, districts.distcode, dffull, radius, Ndist)
    data = (; Y, D, from, to, AP, P, poporig)

    @model function model(Y, from, to, D, AP, P, DM, N, Nfull, fromfull, tofull,
                          Dfull, nfrom, radius, scaleo)
        α ~ Normal(-8, 1)
        β ~ Gamma(1, 1)
        γ_raw ~ Gamma(15, .2); γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1); ϕ = ϕ_raw / 100
        δ_raw ~ Gamma(10, 1); δ = δ_raw / 100

        T = eltype(γ)  # to get dual data type for AD
        denom = zeros(T, nfrom)
        ps = Vector{T}(undef, N)

        if scaleo == 1
            desmat = [exp(desirability(P[i], DM[i, j], ϕ, δ, γ))
                      for i in 1:Ndist, j in 1:Ndist]
        elseif scaleo == "β"
            desmat = [exp(desir(i, j, P[i], DM[i, j], ϕ, δ, γ, β))
                      for i in 1:Ndist, j in 1:Ndist]
        end
        denom = [sum(desmat[i,j] for j in 1:Ndist) for i in 1:Ndist]

        @inbounds for i in 1:N
            ps[i] = AP[from[i]] * exp(α) * desmat[from[i], to[i]] / denom[from[i]]
        end
        Y .~ Poisson.(ps)
        return ps
    end

    mdl = model(Y, from, to, D, AP, P, DM, N, Nfull, fromfull,
                tofull, Dfull, nfrom, radius, scaleo)
    lb = [-20, -100, 10, 0, 1]
    ub = [0, 10, 100, 99, 100]

    return (; mdl, lb, ub, data)
end

desirability(P, D, ϕ, δ, γ) = P + log((ϕ + (1 - ϕ) / (D + δ)^γ))
desir(i, j, P, D, ϕ, δ, γ, β) = i == j ? exp(β) * desirability(P, D, ϕ, δ, γ) : desirability(P, D, ϕ, δ, γ)

function make_distmat(from, to, dffull, radius, Ndist)
    df = DataFrame(from = repeat(from, inner = Ndist),
                   to = repeat(to, outer = Ndist))
    df2 = DataFrame(from = dffull.fromdist, to = dffull.todist,
                    dist = dffull.dist)
    leftjoin!(df, df2, on = [:from, :to])
    df3 = DataFrame(from = from, to = from, radius = radius)
    leftjoin!(df, df3, on = [:from, :to])
    df.dist = coalesce.(df.dist, df.radius)
    return reshape(df.dist, Ndist, Ndist)' ./ 100.0
end

function kon(data::NamedTuple; scaleo)
    df        = data.df
    districts = data.districts
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
                          Dfull, nfrom, radius, scaleo)
        α ~ Normal(-8, 1)
        β ~ Gamma(1, 1)
        γ_raw ~ Gamma(15, .2); γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1); ϕ = ϕ_raw / 100
        δ_raw ~ Gamma(10, 1); δ = δ_raw / 100

        T = eltype(γ)  # to get dual data type for AD
        att   = Vector{T}(undef, N)
        denom = zeros(T, nfrom)
        ps = Vector{T}(undef, N)

        if scaleo == "β"
            c = exp(β)
        elseif scaleo == 1
            c = 1
        end

        @inbounds for i in 1:N
            att[i] = exp(desirability(P[to[i]], D[i], ϕ, δ, γ))
        end

        @inbounds for i in 1:Nfull
            denom[fromfull[i]] += exp(
                desirability(P[tofull[i]], Dfull[i], ϕ, δ, γ))
        end

        @inbounds for i in 1:nfrom
            denom[i] += c * desirability(AP[i], radius[i], ϕ, δ, γ)
        end

        @inbounds for i in 1:N
            ps[i] = AP[from[i]] * exp(α) * att[i] / denom[from[i]]
        end
        Y .~ Poisson.(ps)
        return ps
    end

    mdl = model(Y, from, to, D, AP, P, N, Nfull, fromfull,
                tofull, Dfull, nfrom, radius, scaleo)
    lb = [-20, -100, 10, 0, 1]
    ub = [0, 10, 100, 99, 100]

    return (; mdl, lb, ub, data)
end
