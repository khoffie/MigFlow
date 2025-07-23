getdeviance(r) = deviance(r.mdl.mdl.args.Y, r.prd)[1]

function reorder(results)
    years = [Int(r.chn[:year].data[1]) for r in results]
    return results[sortperm(years)]
end

function loopstruct(s, f, ages = nothing, years = nothing)
    if isnothing(ages); ages = keys(s); end
    if isnothing(years); years = keys(data[ages[1]]); end
    res = Matrix{Any}(undef, length(ages), length(years))
    for i in eachindex(ages)
        a = ages[i]
        for j in eachindex(years)
            y = years[j]
            res[i, j] = f(getfield(getfield(s, a), y))
        end
    end
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
plotcoef(df, c, g, lw = 2) = (Plots.plot(df.year, df[!, c], group = df[!, g], title = c, linewidth = lw))
