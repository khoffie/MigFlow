function gravity(data::NamedTuple; ds = 100, trunc, norm = false)

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
    Ndist      = length(districts.distcode)
    N          = length(Y)
    meta       = MetaData(model = "gravity", age = age, year = year)

    transition(P, D, δ, κ, γ) = κ * P + log(1 / (D + .01) ^ γ

    @model function model(Y::Vector{Int}, from::Vector{Int}, to::Vector{Int},
                          A::Vector{Int}, P::Vector{Float64}, D::Vector{Float64},
                          N::Int, trunc::Bool)

        α_raw ~ Normal(-5, 1);   α = α_raw
        δ_raw ~ Gamma(10, 1);    δ = δ_raw / 10
        κ_raw ~ Gamma(10, 1);    κ = κ_raw / 10
        γ_raw ~ Gamma(15, 0.2);  γ = γ_raw / 10

        T = eltype(γ)
        λ = Vector{T}(undef, N)

        @inbounds for i in 1:N
            λ[i] = A[i]^δ * exp(α + transition(P[to[i]], D[i], δ, κ, γ))
        end
        Y ~ product_distribution(trunc ? TruncatedPoisson.(λ) : Poisson.(λ))
        return trunc ? λ ./ (1 .- exp.(-λ)) : λ
    end

    mdl = model(Y, from, to, A, P, D, N, trunc)
    lb = [-12.0, 1.0, 1.0, 10.0]
    ub = [-0.0 , 40.0, 40.0, 40.0]
    return ModelWrapper(mdl, lb, ub, meta)
end
