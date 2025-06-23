function insensitive(P, α, ϕ)
    return α * P * ϕ
end

function sensitive(P, D, α, γ, δ, ϕ, ds)
    return α * P * ((1 - ϕ) / ((D + δ) ^ γ))
end

function fractioninsensitive(P, D, mles, ds = 100)
    P = P / 153000; D =  D / ds
    α = exp(mles["α_raw"])
    γ = mles["γ_raw"] / 10
    δ = mles["δ_raw"] / 100
    ϕ = mles["ϕ_raw"] / 100
    insens = insensitive.(P, α, ϕ)
    sens = sensitive.(P, D, α, γ, δ, ϕ, ds)
    all = sum(insens + sens)
    frac = sum(insens) / all
    return (; insens, sens, frac, all)
end

function plotinsensitive(insens, sens, dist)
    frac = insens ./ (sens .+ insens)
    plot()
    smoother!(dist, frac)
end

# test = fractioninsensitive(data.df.topop, data.df.dist, out.out, 100)
# plotinsensitive(test.insens, test.sens, data.df.dist)
