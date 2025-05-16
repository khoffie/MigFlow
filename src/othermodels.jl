
function distonly(data::NamedTuple)
    df = data.df
    districts = data.districts
    Y       = df.flows
    from    = levelcode.(categorical(df.fromdist))
    to      = levelcode.(categorical(df.todist))
    AP      = unique(df, :fromdist).frompop ## population of agegroup
    poporig = districts.pop
    P       = log.(districts.pop ./ 153000) # median pop
    D       = df.dist ./ 100.0
    N       = length(Y)
    data = (; Y, D, from, to, AP, P, poporig)

    @model function model(Y, from, to, D, AP, P, N)
        α ~ Normal(-8, 1)
        γ_raw ~ Gamma(15, .2); γ = γ_raw / 10
        ϕ_raw ~ Gamma(10, 1); ϕ = ϕ_raw / 100
        δ_raw ~ Gamma(10, 1); δ = δ_raw / 100

        T = eltype(γ)  # to get dual data type for AD
        ps = Vector{T}(undef, N)

        @inbounds for i in 1:N
            ps[i] = AP[from[i]] * exp(desirability(P[to[i]], D[i], α, ϕ, δ, γ))
        end

        Y .~ Poisson.(ps)

        return ps
    end

    mdl = model(Y, from, to, D, AP, P, N)
    lb = [-20, 10, 0, 1]
    ub = [0, 100, 99, 100]

    return (; mdl, lb, ub, data)
end

desirability(P, D, ϕ, δ, γ) = P + log((ϕ + (1 - ϕ) / (D + δ)^γ))

# function t(x, c, m = 100.0, e = 0.1)
#     B = (c + 1) / ((m + e)^(c+1) - e^(c+1))
#     return B * (x + e)^c
# end

# function gravity(data::NamedTuple)
#     flows = data.flows
#     fp = data.frompop
#     tp = data.topop
#     di = data.dist
#     ds = data.distscale

#     @model function model(flows, fp, tp, di, ds = 100)
#         a ~ Normal(-8, 1)
#         c ~ Gamma(15, .2) ## we want mean ~ 3 but with less variance than Gamma(3, 1)
#         s ~ Gamma(1, 1)
#         t ~ Gamma(1, 1)

#         c = c / 10 # to offest smaller variance of c
#         di = di ./ ds
#         tp = tp ./ 153000 # median pop

#         pop = log.(tp.^t)
#         dist = log.(1 ./ (di).^c)
#         preds = fp.^s .* logistic.(a .+ pop .+ dist)

#         flows ~ arraydist(Poisson.(preds))
#         return preds
#     end

#     mdl = model(flows, fp, tp, di, ds)
#     lb = [-20, 10, 0, 0]
#     ub = [0, 100, 5, 5]
#     return (; mdl, lb, ub)
# end
