function norm(data::NamedTuple; norm, type)
    @assert type ∈ ["joint", "conditional"] "Invalid type $type. Must be 'joint' or 'conditional'."
    @assert norm ∈ ["none", "both", "origin", "destination"] "Invalid norm $norm. Must be 'none', 'both', 'origin' or 'destination'."
    df        = sort(data.df, :fromdist)
    districts = sort(data.districts, :distcode)
    ds        = 100

    Y       = df.flows
    from    = lc(df.fromdist)
    to      = lc(df.todist)
    A       = genfrompop(df, type)
    P       = districts.pop ./ 153000 # median topop
    poporig = districts.pop
    D       = df.dist  ./ ds
    Ndist   = length(districts.distcode)
    N       = length(Y)
    radius  = fradius.(districts.pop, districts.density)
    data    = (; Y, D, from, to, A,  P, poporig)

    @model function model(Y, from, to, A, P, D, Ndist, N, radius, norm)
        α      ~ Normal(-5, 1)
        β      ~ Gamma(1, 1);     ## b  = b_raw / 100
        γ_raw  ~ Gamma(15, 0.2);  γ = γ_raw / 10
        ϕ_raw  ~ Gamma(10, 1.0);  ϕ = ϕ_raw / 100
        δ_raw  ~ Gamma(10, 1.0);  δ = δ_raw / 100

        T = eltype(γ)  # to get dual data type for AD
        att   = Vector{T}(undef, N)
        denom = zeros(T, Ndist)
        ps = Vector{T}(undef, N)

        @inbounds for i in 1:N
            att[i] = desirability(P[to[i]], D[i], γ, δ, ϕ)
        end

        denom = normalize(norm, denom, N, from, att, Ndist, β,
                          desirability, P, radius, γ, δ, ϕ)

        @inbounds for i in 1:N
               ps[i] = A[i] * exp(α) * (att[i] / denom[from[i]])
        end

        Y .~ Poisson.(ps)
        return ps
    end

    mdl = model(Y, from, to, A, P, D, Ndist, N,
                radius, norm)
    lb = [-20.0, -100.0, 10.0, 0.0, 1.0]
    ub = [20.0, 10.0, 100.0, 99.0, 100.0]
    return (; mdl, lb, ub, data)
end

desirability(P, D, γ, δ, ϕ) = P * (ϕ + (1 - ϕ) / ((D + δ) ^ γ))
fradius(P, ρ) = sqrt((P / ρ) / 2π)
lc(x) = levelcode.(categorical(x))

function genfrompop(df, type)
    type == "joint" && return df.frompop
    df2 = combine(groupby(data.df, :fromdist), :flows => sum)
    return leftjoin(df, df2, on = :fromdist).flows_sum
end

function normalize(norm, denom, N, from, att, Ndist, β, desf, P, radius, γ, δ, ϕ)
    norm == "none" && return ones(N)
    if norm in ("destination", "both")
        @inbounds for i in 1:N
            denom[from[i]] += att[i]
        end
    end
    if norm in ("origin", "both")
        @inbounds for i in 1:Ndist
            denom[i] += exp(β) * desf(P[i], radius[i], γ, δ, ϕ)
        end
    end
    return denom
end
