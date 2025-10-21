function baseflownormalized(data::NamedTuple; kdens = 1.5, kgeo = 1.5, ndc = 4, ngcx = 2, ds = 100, trunc)
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
    meta       = MetaData(model = "baseflownormalized", age = age, year = year)

    function desirability(P, D, Q, Gfrom, Gto, γ, ϕ)
        P + log(ϕ + (1 - ϕ) / ((D + 0.01) ^ γ)) + Q + (Gto - Gfrom)
    end

    @model function model(Y::Vector{Int}, from::Vector{Int}, to::Vector{Int},
                          A::Vector{Int}, P::Vector{Float64}, D::Vector{Float64},
                          R::Vector{Float64}, Ndist::Int, N::Int, radius::Vector{Float64},
                          xcoord::Vector{Float64}, ycoord::Vector{Float64},
                          xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64,
                          Rmin::Float64, Rmax::Float64,
                          ndc::Int, ngcx::Int, ngcy::Int,
                          cx, cy, rbf_scale, cxgeo, cygeo, geo_scale, anchor, trunc)

        α_raw ~ Normal(-5, 1);   α = α_raw
        β_raw ~ Gamma(1, 1);     β = β_raw
        γ_raw ~ Gamma(15, 0.2);  γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1.0);  ϕ = ϕ_raw / 100
        # δ_raw ~ Gamma(10, 1.0);  δ = δ_raw / 100
        ζ_raw ~ StMvN(ndc, 10.0);  ζ = coefmat(ζ_raw / 10)
        η_raw ~ StMvN(ngcx * ngcy, 10.0);  η = coefmat(η_raw / 10, ngcx, ngcy)

        T = eltype(γ)  # to get dual data type for AD
        ps = Vector{T}(undef, N)

        des_by_origin = Vector{Vector{T}}(undef, Ndist)
        logden = zeros(T, Ndist)

        for o in 1:Ndist
            des_by_origin[o] = Vector{T}()
        end
        @inbounds for i in 1:N
            push!(des_by_origin[from[i]],
                  desirability(P[to[i]], D[i],
                               interp_anchor(R[from[i]], R[to[i]], ζ, cx, cy, rbf_scale, anchor),
                               interp(xcoord[from[i]], ycoord[from[i]], η, cxgeo, cygeo, geo_scale),
                               interp(xcoord[to[i]], ycoord[to[i]], η, cxgeo, cygeo, geo_scale),
                               γ, ϕ))
        end
        @inbounds for i in 1:Ndist
            push!(des_by_origin[i],
                  desirability(P[i], β * radius[i],
                               interp_anchor(R[i], R[i], ζ, cx, cy, rbf_scale, anchor),
                               0, 0,
                               # interp(xcoord[i], ycoord[i], η, cxgeo, cygeo, geo_scale),
                               # interp(xcoord[i], ycoord[i], η, cxgeo, cygeo, geo_scale),
                               γ, ϕ))
        end

        for o in 1:Ndist
            logden[o] = logsumexp(des_by_origin[o])
        end

        @inbounds for i in 1:N
            ps[i] = A[i] * exp(α +
                desirability(P[to[i]], D[i],
                             interp_anchor(R[from[i]], R[to[i]], ζ, cx, cy, rbf_scale, anchor),
                             interp(xcoord[from[i]], ycoord[from[i]], η, cxgeo, cygeo, geo_scale),
                             interp(xcoord[to[i]], ycoord[to[i]], η, cxgeo, cygeo, geo_scale),
                              γ, ϕ)  - logden[from[i]])
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

    mdl = model(Y, from, to, A, P, D, R, Ndist, N, radius,
                xcoord, ycoord, xmin, xmax, ymin, ymax, Rmin, Rmax,
                ndc, ngcx, ngcy, cx, cy, rbf_scale, cxgeo, cygeo,
                geo_scale, anchor, trunc)
    lb, ub = boundbeta(age, ndc, ngcx, ngcy)
    return ModelWrapper(mdl, lb, ub, meta)
end
