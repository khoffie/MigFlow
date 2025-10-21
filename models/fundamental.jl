function fundamental(data::NamedTuple; ds = 100, trunc, norm)

    df        = sort(data.df, [:fromdist, :todist])
    districts = sort(data.districts, :distcode)
    age       = unique(df.agegroup)[1]
    year      = unique(df.year)[1]

    Y          = df.flows
    from       = lc(df.fromdist)
    to         = lc(df.todist)
    A          = genfrompop(df, "joint")
    # ϵ = 10000, 153000 is median pop in 2017,
    # using actual median(districts.pop) per year leads to failed optimization
    P          = log.((districts.pop .+ 1000) ./ 153000)
    poporig    = districts.pop
    D          = fdist.(df.dist, ds)
    N          = length(Y)
    Ndist      = length(districts.distcode)
    radius     = fradius.(districts.pop, districts.density, ds)
    meta       = MetaData(model = modelname("fundamental", trunc, norm),
                          age = age, year = year)

    transition(P, D, γ, ϕ) = P + log(ϕ + (1 - ϕ) / (D + .01) ^ γ)

    @model function model(Y::Vector{Int}, from::Vector{Int}, to::Vector{Int},
                          A::Vector{Int}, P::Vector{Float64}, D::Vector{Float64},
                          N::Int, Ndist::Int, radius::Vector{Float64}, trunc::Bool, norm)

        α_raw ~ Normal(-5, 1);   α = α_raw
        if norm
            β_raw ~ Gamma(1, 1); β = β_raw
        end
        γ_raw ~ Gamma(15, 0.2);  γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1.0);  ϕ = ϕ_raw / 100

        T = eltype(γ)
        λ = Vector{T}(undef, N)
        Ω = zeros(T, Ndist)

        if norm
            Ω_origin = Vector{Vector{T}}(undef, Ndist)
            for o in 1:Ndist
                Ω_origin[o] = Vector{T}()
            end
            @inbounds for i in 1:N
                push!(Ω_origin[from[i]], transition(P[to[i]], D[i], γ, ϕ))
            end
            @inbounds for i in 1:Ndist
                push!(Ω_origin[i], transition(P[i], β * radius[i], γ, ϕ))
            end
            for o in 1:Ndist
                Ω[o] = logsumexp(Ω_origin[o])
            end
        end

        @inbounds for i in 1:N
            λ[i] = A[i] * exp(α + transition(P[to[i]], D[i], γ, ϕ) - Ω[from[i]])
        end
        Y ~ product_distribution(trunc ? TruncatedPoisson.(λ) : Poisson.(λ))
        return trunc ? λ ./ (1 .- exp.(-λ)) : λ
    end

    mdl = model(Y, from, to, A, P, D, N, Ndist, radius, trunc, norm)
    lb, ub = bound(age, 0, 0, 0, norm)
    return ModelWrapper(mdl, lb, ub, meta)
end
