function decaytime(y₀, yₜ, r, steptime = 5)
    λ = log(1 - r)
    s = log(yₜ / y₀) / λ
    t = s * steptime / 60
    out = (steps = round(s, digits = 2), time = round(t, digits = 2))
    return out
end

function chainnames(path, temp = nothing)
    fchains = [f for f in readdir(path) if contains(f, "germchain")]
    if temp == nothing
        fchains = fchains[.!contains.(fchains, ".0.csv")]
    elseif temp == "all"
        fchains
    else
        fchains = fchains[contains.(fchains, "$(temp).0.csv")]
    end
    return fchains
end
