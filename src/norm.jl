function norm(data::NamedTuple; kdens = 1.5, kgeo = 1.5, ndc = 4, ngcx = 2, normalize = true, ds = 100)

    df        = sort(data.df, [:fromdist, :todist])
    districts = sort(data.districts, :distcode)
    dffull    = sort(data.dffull, [:fromdist, :todist])
    age       = unique(df.agegroup)[1]
    year      = unique(df.year)[1]

    Y          = df.flows
    from       = lc(df.fromdist)
    to         = lc(df.todist)
    A          = genfrompop(df, "joint")
    P          = log.(districts.pop ./ 153000) # median topop
    poporig    = districts.pop
    D          = fdist.(df.dist, ds)
    R          = scale_to_unit(log.(districts.density ./ median(districts.density)))
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
    ## Germany's height is about 1.33 times it's with. Of the width,
    ## we only use 50% as space for the grid, of the height we use
    ## 70%. For equidistant points we thus need 1.86 = 1.33 * .7 / .5 as many
    ## different y points as x points
    ngcy = Int(round(ngcx * 1.86, digits = 0))
    cxgeo      = [range(-0.6, 0.4, ngcx);]
    cygeo      = [range(-0.9, 0.6, ngcy);]
    geo_scale  = rbfscale(cxgeo, cygeo, kgeo)
    fromfull   = lc(dffull.fromdist)
    tofull     = lc(dffull.todist)
    Dfull      = fdist(dffull.dist, ds)
    Nfull      = length(fromfull)
    data       = (; Y, D, from, to, A,  P, R, districts.distcode, poporig,
                  ndc, ngcx, ngcy, Rmin, Rmax, xcoord, ycoord, age, year,
                  cx, cy, cxgeo, cygeo, kdens, kgeo)

    @model function model(Y::Vector{Int}, from::Vector{Int}, to::Vector{Int},
                          A::Vector{Int}, P::Vector{Float64}, D::Vector{Float64},
                          R::Vector{Float64}, Ndist::Int, N::Int, radius::Vector{Float64},
                          xcoord::Vector{Float64}, ycoord::Vector{Float64},
                          xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64,
                          Rmin::Float64, Rmax::Float64,
                          fromfull::Vector{Int}, tofull::Vector{Int},
                          Dfull::Vector{Float64}, Nfull::Int,
                          ndc::Int, ngcx::Int, ngcy::Int,
                          normalize::Bool, cx, cy, rbf_scale, cxgeo, cygeo, geo_scale)

        α_raw ~ Normal(-5, 1);   α = α_raw
        β_raw ~ Gamma(1, 1);     β = β_raw
        γ_raw ~ Gamma(15, 0.2);  γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1.0);  ϕ = ϕ_raw / 100
        δ_raw ~ Gamma(10, 1.0);  δ = δ_raw / 100
        ζ_raw ~ StMvN(ndc, 10.0);  ζ = coefmat(ζ_raw / 10)
        η_raw ~ StMvN(ngcx * ngcy, 10.0);  η = coefmat(η_raw / 10, ngcx, ngcy)

        T = eltype(γ)  # to get dual data type for AD
        denom = zeros(T, Ndist)
        ps = Vector{T}(undef, N)

        if normalize
            @inbounds for i in 1:Nfull
                denom[fromfull[i]] += desirability(P[tofull[i]], Dfull[i],
                                                   interpolant(rbf, R[fromfull[i]], R[tofull[i]], ζ, cx, cy, rbf_scale),
                                                   interpolant(rbf, xcoord[fromfull[i]], ycoord[fromfull[i]], η, cxgeo, cygeo, geo_scale),
                                                   interpolant(rbf, xcoord[tofull[i]], ycoord[tofull[i]], η, cxgeo, cygeo, geo_scale),
                                                   γ, δ, ϕ)
            end
            @inbounds for i in 1:Ndist
                denom[i] += desirability(P[i], β * radius[i],
                                         interpolant(rbf, R[i], R[i], ζ, cx, cy, rbf_scale),
                                         1, 1, γ, δ, ϕ)
            end
        else
            fill!(denom, one(T))
        end

        @inbounds for i in 1:N
            ps[i] = A[i] * exp(α +
                desirability(P[to[i]], D[i],
                             interpolant(rbf, R[from[i]], R[to[i]], ζ, cx, cy, rbf_scale),
                             interpolant(rbf, xcoord[from[i]], ycoord[from[i]], η, cxgeo, cygeo, geo_scale),
                             interpolant(rbf, xcoord[to[i]], ycoord[to[i]], η, cxgeo, cygeo, geo_scale),
                              γ, δ, ϕ) / denom[from[i]])
        end
        Y ~ product_distribution(Poisson.(ps))
        return ps
    end

    mdl = model(Y, from, to, A, P, D, R, Ndist, N, radius,
                xcoord, ycoord, xmin, xmax, ymin, ymax, Rmin, Rmax,
                fromfull, tofull, Dfull, Nfull, ndc, ngcx, ngcy, normalize,
                cx, cy, rbf_scale, cxgeo, cygeo, geo_scale)
    lb = [-20.0, -100.0, 10.0, 1.0, 1.0, fill(-100, ndc)..., fill(-100, ngcx * ngcy)...]
    ub = [20.0, 100.0, 100.0, 99.0, 99.0, fill(100, ndc)..., fill(100, ngcx * ngcy)...]
    return (; mdl, lb, ub, data)
end

# function desirability(P, D, Q, Gfrom, Gto, γ, δ, ϕ)
#     P * (ϕ + (1 - ϕ) / ((D + δ) ^ γ) * Q * (Gto / Gfrom))
# end

function desirability(P, D, Q, Gfrom, Gto, γ, δ, ϕ)
    P + log(ϕ + (1 - ϕ) / ((D + δ) ^ γ)) + Q + (Gto - Gfrom)
end

fdist(D, ds) = D / ds
fradius(P, ρ, ds) = sqrt((P / ρ) / 2π) / ds
lc(x) = levelcode.(categorical(x))
StMvN(n, σ) = MvNormal(zeros(n), fill(σ, n))
function coefmat(c)
    s = Int(sqrt(length(c)))
    return reshape(c, (s, s))
end
function coefmat(c, nx::Int, ny::Int)
    return reshape(c, (nx, ny))
end

function genfrompop(df, type)
    type == "joint" && return df.frompop
    df2 = combine(groupby(data.df, :fromdist), :flows => sum)
    return leftjoin(df, df2, on = :fromdist).flows_sum
end
