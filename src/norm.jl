function norm(data::NamedTuple; ndc = 1, ngc = 1, normalize = true)

    df        = sort(data.df, [:fromdist, :todist])
    districts = sort(data.districts, :distcode)
    dffull    = sort(data.dffull, [:fromdist, :todist])
    ds        = 100

    Y          = df.flows
    from       = lc(df.fromdist)
    to         = lc(df.todist)
    A          = genfrompop(df, "joint")
    P          = districts.pop ./ 153000 # median topop
    poporig    = districts.pop
    D          = fdist.(df.dist, ds)
    R          = log.(districts.density ./ median(districts.density))
    Rmin, Rmax = extrema(R)
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
                  ndc, ngc, Rmin, Rmax, xcoord, ycoord)

    @model function model(Y, from, to, A, P, D, R, Ndist, N, radius,
                          xcoord, ycoord, xmin, xmax, ymin, ymax, Rmin, Rmax,
                          fromfull, tofull, Dfull, Nfull, ndc, ngc,
                          normalize)

        α_raw  ~ Normal(-5, 1);   α = exp(α_raw)
        β_raw  ~ Gamma(1, 1);     β = exp(β_raw)
        γ_raw  ~ Gamma(15, 0.2);  γ = γ_raw / 10
        ϕ_raw  ~ Gamma(10, 1.0);  ϕ = ϕ_raw / 100
        δ_raw  ~ Gamma(10, 1.0);  δ = δ_raw / 100
        ζ_raw ~ StMvN(ndc, 10); ζ = ζ_raw / 100; ζ[1] = 0.0 # cheby intercept, ensure heatmap is not elevated
        η_raw ~ StMvN(ngc, 10); η = η_raw / 100; η[1] = 0.0

        T = eltype(γ)  # to get dual data type for AD
        denom = zeros(T, Ndist)
        ps = Vector{T}(undef, N)

        Q = log.( 1 .+ exp.(defdensitycheby(ζ, Rmin, Rmax).(R[from], R[to])))
        G = exp.(defgeocheby(η, xmin, xmax, ymin, ymax).(xcoord, ycoord))

        if normalize
            Qfull = exp.(defdensitycheby(ζ, Rmin, Rmax).(R[fromfull], R[tofull]))
            @inbounds for i in 1:Nfull
                denom[fromfull[i]] += desirability(P[tofull[i]], Dfull[i], Qfull[i],
                                                   G[fromfull[i]], G[tofull[i]], γ, δ, ϕ)
            end
            @inbounds for i in 1:Ndist
                denom[i] += desirability(P[i], β * radius[i],
                                         exp.(defdensitycheby(ζ, Rmin, Rmax).(R[i], R[i])),
                                         1, 1, γ, δ, ϕ)
            end
        else
            denom = ones(T, Ndist)
        end

        @inbounds for i in 1:N
            ps[i] = A[i] * α *
                (desirability(P[to[i]], D[i], Q[i],
                              G[from[i]], G[to[i]], γ, δ, ϕ) / denom[from[i]])
        end
        Y .~ Poisson.(ps)
        return ps
    end

    mdl = model(Y, from, to, A, P, D, R, Ndist, N, radius,
                xcoord, ycoord, xmin, xmax, ymin, ymax, Rmin, Rmax,
                fromfull, tofull, Dfull, Nfull, ndc, ngc, normalize)
    lb = [-20.0, -100.0, 10.0, 0.0, 1.0, fill(-100, ndc)..., fill(-100, ngc)...]
    ub = [20.0, 100.0, 100.0, 99.0, 100.0, fill(100, ndc)..., fill(100, ngc)...]
    return (; mdl, lb, ub, data)
end

## desirability(P, D, γ, δ, ϕ) = P * (ϕ + (1 - ϕ) / ((D + δ) ^ γ))
function desirability(P, D, Q, Gfrom, Gto, γ, δ, ϕ)
    P * (ϕ + (1 - ϕ) / ((D + δ) ^ γ) * Q * exp(Gto - Gfrom))
end

fdist(D, ds) = D / ds
fradius(P, ρ, ds) = sqrt((P / ρ) / 2π) / ds
lc(x) = levelcode.(categorical(x))
StMvN(n, σ) = MvNormal(zeros(n), fill(σ, n))

function genfrompop(df, type)
    type == "joint" && return df.frompop
    df2 = combine(groupby(data.df, :fromdist), :flows => sum)
    return leftjoin(df, df2, on = :fromdist).flows_sum
end
