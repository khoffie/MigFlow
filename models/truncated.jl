function truncated(data::NamedTuple; ds = 100)

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
    meta       = MetaData(model = "truncated", age = age, year = year)

    @model function model(Y::Vector{Int}, from::Vector{Int}, to::Vector{Int},
                          A::Vector{Int}, P::Vector{Float64}, D::Vector{Float64},
                          N::Int)

        α_raw ~ Normal(-5, 1);   α = α_raw
        γ_raw ~ Gamma(15, 0.2);  γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1.0);  ϕ = ϕ_raw / 100

        T = eltype(γ)  # to get dual data type for AD
        ps = Vector{T}(undef, N)

        @inbounds for i in 1:N
            ps[i] = A[i] * exp(α + P[to[i]] + log(ϕ + (1 - ϕ) / (D[i] + .01) ^ γ))
        end
        Y ~ product_distribution(Poisson.(ps))
        return ps
    end

    mdl = model(Y, from, to, A, P, D, N)
    lb = [-12.0, 10.0, 1.0]
    ub = [-5.0 , 40.0, 50.0]
    return ModelWrapper(mdl, lb, ub, meta)
end
