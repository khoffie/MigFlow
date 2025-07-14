getdeviance(r) = deviance(r.mdl.mdl.args.Y, returned(r.mdl.mdl, r.chn)[1])

function reorder(results)
    years = [Int(r.chn[:year].data[1]) for r in results]
    return results[sortperm(years)]
end

function loopstruct(s, f, ages = nothing, years = nothing)
    if isnothing(ages)
        ages = fieldnames(typeof(s))
    end
    if isnothing(years)
        years = fieldnames(typeof(getfield(s, ages[1])))
    end
    res = [f(getfield(getfield(s, a), y)) for a in ages, y in years]
    return res
end

function extract_params(result)
    params = result.chn.name_map.parameters
    df = DataFrame([Symbol(p) => result.chn[p].data[1] for p in params])
    df.deviance .= getdeviance(result)
    df.year .= result.chn[:year].data[1]
    df.age .= recodeage(Int(result.chn[:age].data[1]))
    first = ["age", "year", "α_raw", "γ_raw", "ϕ_raw", "lp"]
    last = setdiff(names(df), first)
    select!(df, vcat(first, last))
    return df
end

plotcoef(df, c) = (plot(df.group, df[!, c], title = c); scatter!(df.group, df[!, c]))
plotcoef(df, c, g) = (plot(df.year, df[!, c], group = df[!, g], title = c))
