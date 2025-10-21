function normalized(data::NamedTuple; ds = 100, trunc = false)
    ## this normalization here will only work with full data, so use
    ## p=1.0 if generating the data
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
    D          = fdist.(df.dist, ds)
    N          = length(Y)
    Ndist      = length(districts.distcode)
    radius     = fradius.(districts.pop, districts.density, ds)
    meta       = MetaData(model = "normalized", age = age, year = year)

    desirability(P, D, γ, ϕ) = P + log(ϕ + (1 - ϕ) / ((D + 0.01) ^ γ))

    @model function model(Y::Vector{Int}, from::Vector{Int}, to::Vector{Int},
                          A::Vector{Int}, P::Vector{Float64}, D::Vector{Float64},
                          Ndist::Int, N::Int, radius::Vector{Float64}, trunc)

        α_raw ~ Normal(-5, 1);   α = α_raw
        β_raw ~ Gamma(1, 1);     β = β_raw
        γ_raw ~ Gamma(15, 0.2);  γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1.0);  ϕ = ϕ_raw / 100

        T = eltype(γ)  # to get dual data type for AD
        denom = zeros(T, Ndist)
        ps = Vector{T}(undef, N)
        des_by_origin = Vector{Vector{T}}(undef, Ndist)
        logden = zeros(T, Ndist)

        for o in 1:Ndist
            des_by_origin[o] = Vector{T}()
        end
        @inbounds for i in 1:N
            push!(des_by_origin[from[i]], desirability(P[to[i]], D[i], γ, ϕ))
        end
        @inbounds for i in 1:Ndist
            push!(des_by_origin[i], desirability(P[i], β * radius[i], γ, ϕ))
        end

        for o in 1:Ndist
            logden[o] = logsumexp(des_by_origin[o])
        end

        @inbounds for i in 1:N
            ps[i] = A[i] * exp(α + desirability(P[to[i]], D[i], γ, ϕ) - logden[from[i]])
        end
        if !trunc
            Y ~ product_distribution(Poisson.(ps))
            preds = ps
        elseif trunc
            Y ~ product_distribution(TruncatedPoisson.(ps))
            preds = ps ./ (1 .- exp.(-ps))
        end
        return preds
    end

    mdl = model(Y, from, to, A, P, D, Ndist, N, radius, trunc)
    lb = [-12.0, 0.0, 10.0, 1.0]
    ub = [0.0,   100.0, 40.0, 50.0]
    return ModelWrapper(mdl, lb, ub, meta)
end
