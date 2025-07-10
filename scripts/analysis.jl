using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings, Loess
using ADTypes, KernelDensity, Serialization, DynamicPPL, LinearAlgebra
using BenchmarkTools, IterTools, StatProfilerHTML, ReverseDiff
using Suppressor
## using Enzyme

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/diagchain.jl")
include("../src/diagrbf.jl")
include("../src/norm.jl")
include("../src/rbf.jl")

getdeviance(r) = deviance(r.mdl.mdl.args.Y, returned(r.mdl.mdl, r.chn)[1])
function reorder(results)
    years = [Int(r.chn[:year].data[1]) for r in results]
    return results[sortperm(years)]
end

function loopstruct(s, f)
    ages = fieldnames(typeof(s))
    years = fieldnames(typeof(getfield(s, ages[1])))
    res = [@suppress f(getfield(getfield(data, age), year))
        for age in ages, year in years]
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
plotcoef(df, c, g) = (plot(df.group, df[!, c], group = df[!, g], title = c))

function resplots(results, coefdf)
    pgamma = plotcoef(coefdf, :γ_raw)
    palpha = plotcoef(coefdf, :α_raw)
    pphi = plotcoef(coefdf, :ϕ_raw)
    pdeviance = plotcoef(coefdf, :deviance)
    pdens = [r.plts[5] for r in results]
    pgeo = [r.plts[6] for r in results]
    return [palpha, pgamma, pphi, pdeviance, pdens, pgeo]
end

years = vcat(2000:2002, 2004:2017)
files = readdir("./output"; join = true)
ages = ["18to25", "25to30", "30to50", "50to65", "above65", "below18"]

@eval begin
    struct YearResults
        $(Symbol.("y" .* string.(years))...)
    end
end

@eval begin
    struct AgeResults
        $(Symbol.("age" .* ages)...)
    end
end

data = [YearResults(reorder(deserialize(f))...) for f in files];
data = AgeResults(data...);

function loopstruct(s, f, ages = nothing, years = nothing)
    if isnothing(ages)
        ages = fieldnames(typeof(s))
    end
    if isnothing(years)
        years = fieldnames(typeof(getfield(s, ages[1])))
    end
    res = [@suppress f(getfield(getfield(s, age), year))
        for age in ages, year in years]
    return res
end



df = loopstruct(data, extract_params)
df = reduce(vcat, df)



## accidentaly coded above65 as 5 not 6
dfs[5].age .= 6.0
dfs = reduce(vcat, dfs)
dfs.age = recodeage.(Int.(dfs.age))

plot(plotcoef(dfs, :α_raw, :age),
plotcoef(dfs, :γ_raw, :age),
plotcoef(dfs, :ϕ_raw, :age),
plotcoef(dfs, :deviance, :age))

out = [analyze(r.chn, r.mdl) for r in results];
df = coefdf(results)
plts = resplots(out, df)


idx = [1, 2, 3, 15, 16, 17]
plot(plts[1:4]...)
plot(plts[4][idx]..., size = (1200, 900))
plot(plts[5][idx]..., size = (1200, 900))




res = [coefdf(reorder(deserialize(f))) for f in files];



keys(results[1])

typeof(dfs.age)
