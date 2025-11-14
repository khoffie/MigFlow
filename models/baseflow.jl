function baseflow(data::NamedTuple; kdens = 1.5, kgeo = 1.5, ndc = 4, ngcx = 2, ds = 100,
                            trunc, norm)
    ## why sort?
    ## Ndist and radius not needed anymore
    ## poporig also not needed?
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
    R          = scale_to_unit(log.(districts.density ./ median(districts.density)))
    anchor     = median(R)
    Rmin, Rmax = extrema(R)
    cx         = [range(Rmin, Rmax, Int(sqrt(ndc)));]
    cy         = [range(Rmin, Rmax, Int(sqrt(ndc)));]
    rbf_scale  = rbfscale(cx, cy, kdens)
    Ndist      = length(districts.distcode)
    N          = length(Y)
    radius     = fradius.(districts.pop, districts.density, ds)
    xcoord     = scale_to_unit(districts.xcoord)
    ycoord     = scale_to_unit(districts.ycoord)
    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    xlim       = (- 1.0, 1.0)
    ylim       = (- 1.0, 1.0)
    r          = geo_ratio(districts.xcoord, districts.ycoord, xlim, ylim)
    ngcy       = Int(round(r * ngcx, digits = 0))
    cxgeo      = [range(xlim[1], xlim[2], ngcx);]
    cygeo      = [range(ylim[1], ylim[2], ngcy);]
    geo_scale  = rbfscale(cxgeo, cygeo, kgeo)
    meta       = MetaData(model = modelname("baseflow", trunc, norm), age = age, year = year)

    function transition(P, D, Q, Gfrom, Gto, γ, ϕ)
        P + log(ϕ + (1 - ϕ) / ((D + .01) ^ γ)) + Q + (Gto - Gfrom)
    end

    @model function model(Y::Vector{Int}, from::Vector{Int}, to::Vector{Int},
                          A::Vector{Int}, P::Vector{Float64}, D::Vector{Float64},
                          R::Vector{Float64}, Ndist::Int, N::Int, radius::Vector{Float64},
                          xcoord::Vector{Float64}, ycoord::Vector{Float64},
                          xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64,
                          Rmin::Float64, Rmax::Float64,
                          ndc::Int, ngcx::Int, ngcy::Int,
                          cx, cy, rbf_scale, cxgeo, cygeo, geo_scale, anchor, trunc, norm)

        α_raw ~ Normal(-5, 1);   α = α_raw
        if norm
            β_raw ~ Gamma(1, 1); β = β_raw
        end
        γ_raw ~ Gamma(15, 3 / 14);  γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 10 / 9);  ϕ = ϕ_raw / 100
        ζ_raw ~ StMvN(ndc, 10.0);  ζ = coefmat(ζ_raw / 10)
        η_raw ~ StMvN(ngcx * ngcy, 10.0);  η = coefmat(η_raw / 10, ngcx, ngcy)

        T = eltype(γ)
        λ = Vector{T}(undef, N)
        Ω = zeros(T, Ndist)

        if norm
            Ω_origin = Vector{Vector{T}}(undef, Ndist)
            for o in 1:Ndist
                Ω_origin[o] = Vector{T}()
            end
            @inbounds for i in 1:N
                push!(Ω_origin[from[i]],
                      transition(P[to[i]], D[i],
                                 interp_anchor(R[from[i]], R[to[i]], ζ, cx, cy, rbf_scale, anchor),
                                 interp(xcoord[from[i]], ycoord[from[i]], η, cxgeo, cygeo, geo_scale),
                                 interp(xcoord[to[i]], ycoord[to[i]], η, cxgeo, cygeo, geo_scale),
                                 γ, ϕ))
            end
            @inbounds for i in 1:Ndist
                push!(Ω_origin[i],
                      transition(P[i], β * radius[i],
                                 interp_anchor(R[i], R[i], ζ, cx, cy, rbf_scale, anchor),
                                 0, 0, γ, ϕ))
            end
            for o in 1:Ndist
                Ω[o] = logsumexp(Ω_origin[o])
            end
        end

        @inbounds for i in 1:N
            λ[i] = A[i] * exp(α +
                transition(P[to[i]], D[i],
                           interp_anchor(R[from[i]], R[to[i]], ζ, cx, cy, rbf_scale, anchor),
                           interp(xcoord[from[i]], ycoord[from[i]], η, cxgeo, cygeo, geo_scale),
                           interp(xcoord[to[i]], ycoord[to[i]], η, cxgeo, cygeo, geo_scale),
                           γ, ϕ)  - Ω[from[i]])
        end
        Y ~ product_distribution(trunc ? TruncatedPoisson.(λ) : Poisson.(λ))
        return trunc ? λ ./ (1 .- exp.(-λ)) : λ
    end

    mdl = model(Y, from, to, A, P, D, R, Ndist, N, radius,
                xcoord, ycoord, xmin, xmax, ymin, ymax, Rmin, Rmax,
                ndc, ngcx, ngcy, cx, cy, rbf_scale, cxgeo, cygeo,
                geo_scale, anchor, trunc, norm)
    lb, ub = bound(age, ndc, ngcx, ngcy, norm)
    return ModelWrapper(mdl, lb, ub, meta)
end
