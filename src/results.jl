function readresults(patterns::Vector{<:AbstractString}, path="./output")
    results = NamedTuple()
    for pat in patterns
        res, ages = readresults(pat, path)
        m = res[1][1].mdl.meta.model  # extract model name string
        results = merge(results, (; Symbol(m) => res))
    end
    ages = intersect([keys(r) for r in results])[1]
    return results, ages
end

function readresults(pattern::AbstractString, path = "./output")
    years = Symbol.("y" .* string.(vcat(2000:2002, 2004:2017)))
    files = readdir(path; join = true)
    files = files[contains.(files, pattern)]
    ## extracting age group
    ages = string.([s[end] for s in split.(files, "/")])
    ages = "age" .* ([s[2] for s in split.(ages, pattern)])
    ages = Symbol.(replace.(ages, "-" => "to"))

    data = (; (
        age => (; zip(years, reorder(deserialize(file)))...)
        for (age, file) in zip(ages, files)
            )...)
    return data, ages
end

function reorder(results)
    years = [r.mdl.meta.year for r in results]
    return results[sortperm(years)]
end

function loopstruct(s, f, models = nothing, ages = nothing, years = nothing)
    if isnothing(models); models = keys(s); end
    if isnothing(ages); ages = keys(s[models[1]]); end
    if isnothing(years); years = keys(s[models[1]][ages[1]]); end
    res = Array{Any}(undef, length(models), length(ages), length(years))
    for k in eachindex(models)
        m = models[k]
        for i in eachindex(ages)
            a = ages[i]
            for j in eachindex(years)
                y = years[j]
                res[k, i, j] = f(getfield(getfield(getfield(s, m), a), y))
            end
        end
    end
    return res
end

function extract_params(r::EstimationResult)
    params = r.chn.name_map.parameters
    df = DataFrame([Symbol(p) => r.chn[p].data[1] for p in params])
    df.deviance .= getdeviance(r)
    df.year .= r.mdl.data.year
    df.agegroup .= r.mdl.data.age
    df.lp .= r.mdl.data.lp
    first = ["year", "agegroup", "lp"]
    last = setdiff(names(df), first)
    select!(df, vcat(first, last))
    return df
end

plotcoef(df, c) = (plot(df.group, df[!, c], title = c); scatter!(df.group, df[!, c]))
plotcoef(df, c, g, lw = 2) = (Plots.plot(df.year, df[!, c], group = df[!, g], title = c, linewidth = lw))
getdeviance(r) = deviance2(r.mdl.mdl.args.Y, r.prd)[1]
