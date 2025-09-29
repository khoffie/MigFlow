function readresults(path = "./output")
    years = Symbol.("y" .* string.(vcat(2000:2002, 2004:2017)))
    files = readdir(path; join = true)
    ## extracting age group
    ages = string.([s[end] for s in split.(files, "/")])
    ages = "age" .* ([s[2] for s in split.(ages, "optim")])
    ages = Symbol.(replace.(ages, "-" => "to"))

    data = (; (
        age => (; zip(years, reorder(deserialize(file)))...)
        for (age, file) in zip(ages, files)
            )...)
    return data, ages
end

function reorder(results)
    years = [r.mdl.data.year for r in results]
    return results[sortperm(years)]
end

function loopstruct(s, f, ages = nothing, years = nothing)
    if isnothing(ages); ages = keys(s); end
    if isnothing(years); years = keys(s[ages[1]]); end
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

function extract_params(r::EstimationResult)
    params = r.chn.name_map.parameters
    df = DataFrame([Symbol(p) => r.chn[p].data[1] for p in params])
    df.deviance .= getdeviance(r)
    df.year .= r.mdl.data.year
    df.age .= r.mdl.data.age
    df.lp .= r.mdl.data.lp
    first = ["year", "age", "α_raw", "γ_raw", "ϕ_raw", "lp"]
    last = setdiff(names(df), first)
    select!(df, vcat(first, last))
    return df
end

plotcoef(df, c) = (plot(df.group, df[!, c], title = c); scatter!(df.group, df[!, c]))
plotcoef(df, c, g, lw = 2) = (Plots.plot(df.year, df[!, c], group = df[!, g], title = c, linewidth = lw))
getdeviance(r) = deviance2(r.mdl.mdl.args.Y, r.prd)[1]
