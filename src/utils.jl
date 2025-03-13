function sample_rows(df::DataFrame, p::AbstractFloat)
    nrows = nrow(df)
    n = Int(floor(nrows * p))
    df = df[Random.shuffle(1:nrows)[1 : n], : ]
    return(df)
end

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

function get_params2(mdl)
    params = DynamicPPL.syms(DynamicPPL.VarInfo(mdl))
    params = collect(string.(params))
    if "kd" in params
        all_kd = ["kd[$i]" for i in 1:36]
        params = vcat(filter(x -> x != "kd", params), all_kd)
    end
    if "desirecoefs" in params
        all_des = ["desirecoefs[$i]" for i in 1:36]
        params = vcat(filter(x -> x != "desirecoefs", params), all_des)
    end

    return params
end
