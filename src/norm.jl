function norm(data::NamedTuple; densscale = 2.0, ndc = 1, ngc = 1, normalize = true, ds = 100)

    df        = sort(data.df, [:fromdist, :todist])
    districts = sort(data.districts, :distcode)
    dffull    = sort(data.dffull, [:fromdist, :todist])
    age       = unique(df.agegroup)[1]
    year      = unique(df.year)[1]

    Y          = df.flows
    from       = lc(df.fromdist)
    to         = lc(df.todist)
    A          = genfrompop(df, "joint")
    P          = districts.pop ./ 153000 # median topop
    poporig    = districts.pop
    D          = fdist.(df.dist, ds)
    R          = scale_to_unit(log.(districts.density ./ median(districts.density)))
    Rmin, Rmax = extrema(R)
    cx         = [range(Rmin, Rmax, Int(sqrt(ndc)));]
    cy         = [range(Rmin, Rmax, Int(sqrt(ndc)));]
    rbf_scale   = rbfscale(cx, cy, densscale)
    Ndist      = length(districts.distcode)
    N          = length(Y)
    radius     = fradius.(districts.pop, districts.density, ds)
    xcoord     = districts.xcoord
    ycoord     = districts.ycoord
    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    fromfull   = lc(dffull.fromdist)
    tofull     = lc(dffull.todist)
    Dfull      = fdist(dffull.dist, ds)
    Nfull      = length(fromfull)
    data       = (; Y, D, from, to, A,  P, districts.distcode, poporig,
                  ndc, ngc, Rmin, Rmax, xcoord, ycoord, age, year)

    @model function model(Y::Vector{Int}, from::Vector{Int}, to::Vector{Int},
                          A::Vector{Int}, P::Vector{Float64}, D::Vector{Float64},
                          R::Vector{Float64}, Ndist::Int, N::Int, radius::Vector{Float64},
                          xcoord::Vector{Float64}, ycoord::Vector{Float64},
                          xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64,
                          Rmin::Float64, Rmax::Float64,
                          fromfull::Vector{Int}, tofull::Vector{Int},
                          Dfull::Vector{Float64}, Nfull::Int, ndc::Int, ngc::Int,
                          normalize::Bool, cx, cy, rbf_scale)

        α_raw  ~ Normal(-5, 1);   α = exp(α_raw)
        β_raw  ~ Gamma(1, 1);     β = β_raw
        γ_raw  ~ Gamma(15, 0.2);  γ = γ_raw / 10
        ϕ_raw  ~ Gamma(10, 1.0);  ϕ = ϕ_raw / 100
        δ_raw  ~ Gamma(10, 1.0);  δ = δ_raw / 100
        ζ_raw  ~ StMvN(ndc, 10);  ζ = ζ_raw / 100
        η_raw ~ StMvN(ngc, 10);   η = η_raw / 100

        T = eltype(γ)  # to get dual data type for AD
        denom = zeros(T, Ndist)
        ps = Vector{T}(undef, N)

        G = exp.(defgeocheby(η, xmin, xmax, ymin, ymax).(xcoord, ycoord))

        if normalize
            @inbounds for i in 1:Nfull
                denom[fromfull[i]] += desirability(P[tofull[i]], Dfull[i],
                                                   interpolant(rbf, R[fromfull[i]], R[tofull[i]], ζ, cx, cy, rbf_scale),
                                                   ## dc(R[fromfull[i]], R[tofull[i]]),
                                                   G[fromfull[i]], G[tofull[i]],
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
            ps[i] = A[i] * α *
                desirability(P[to[i]], D[i],
                             exp(interpolant(rbf, R[from[i]], R[to[i]], ζ, cx, cy, rbf_scale)),
                             G[from[i]], G[to[i]],
                             γ, δ, ϕ) / denom[from[i]]
        end
        Y ~ product_distribution(Poisson.(ps))
        return ps
    end

    mdl = model(Y, from, to, A, P, D, R, Ndist, N, radius,
                xcoord, ycoord, xmin, xmax, ymin, ymax, Rmin, Rmax,
                fromfull, tofull, Dfull, Nfull, ndc, ngc, normalize,
                cx, cy, rbf_scale)
    lb = [-20.0, -100.0, 10.0, 1.0, 1.0, fill(-200, ndc)..., fill(-100, ngc)...]
    ub = [20.0, 100.0, 100.0, 99.0, 99.0, fill(200, ndc)..., fill(100, ngc)...]
    return (; mdl, lb, ub, data)
end

function desirability(P, D, Q, Gfrom, Gto, γ, δ, ϕ)
    P * (ϕ + (1 - ϕ) / ((D + δ) ^ γ) * Q * (Gto / Gfrom))
end

# function desirability(P, D, Q, Gfrom, Gto, γ, δ, ϕ)
#     P + log(ϕ + (1 - ϕ) / ((D + δ) ^ γ)) + Q + (Gto - Gfrom)
# end


fdist(D, ds) = D / ds
fradius(P, ρ, ds) = sqrt((P / ρ) / 2π) / ds
lc(x) = levelcode.(categorical(x))
StMvN(n, σ) = MvNormal(zeros(n), fill(σ, n))

function genfrompop(df, type)
    type == "joint" && return df.frompop
    df2 = combine(groupby(data.df, :fromdist), :flows => sum)
    return leftjoin(df, df2, on = :fromdist).flows_sum
end
