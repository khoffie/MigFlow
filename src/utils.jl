function decaytime(y₀, yₜ, r, steptime = 5)
    λ = log(1 - r)
    s = log(yₜ / y₀) / λ
    t = s * steptime / 60
    out = (steps = round(s, digits = 2), time = round(t, digits = 2))
    return out
end
