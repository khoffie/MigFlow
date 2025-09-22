age(df, age::Vector{AbstractString}) = filter(:agegroup => n -> n ∈ age, df)
age(df, age::AbstractString) = filter(:agegroup => n -> n == age, df)
year(df, y::Vector{Int64}) = filter(:year => n -> n ∈ y, df)
year(df, y::Int64) = filter(:year => n -> n == y, df)
year(df, xmin::Float64, xmax::Float64) = df[df.year .<= xmax .&& df.year .>= xmin, :]
origin(df, o::Vector{Int64}) = filter(:fromdist => n -> n ∈ o, df)
origin(df, o::Int64) = filter(:fromdist => n -> n == o, df)
destination(df, d) = filter(:todist => n -> n ∈ d, df)
code(df, c) = filter(:distcode => n -> n ∈ c, df)
rmreg(df::DataFrame, col) = filter(col => n -> n != 3159, df) ## data issues
rmreg(df::GeoTable, col) = GeoTable(filter(col => n -> n != 3159, DataFrame(df))) ## data issues

function topn(df, group, col, N = 10)
    dfg = groupby(sort(df, col, rev = true), group)
    dfg = [dfg[i][1:N, [group..., col]] for i in 1 : length(dfg)]
    return reduce(vcat, dfg)
end

function outflux(df, f, group = nothing, name = nothing)
    group == nothing ? g = :fromdist : g = [:fromdist, group...]
    out = combine(DataFrames.groupby(df, g), :flows => f)
    if isnothing(name)
        rename!(out, "flows" * "_$(string(f))" => :flows)
    else
        rename!(out, "flows" * "_$(string(f))" => name)
    end
    return out
end

function influx(df, f, group = nothing)
    group == nothing ? g = :todist : g = [:todist, group...]
    out = combine(DataFrames.groupby(df, g), :flows => f)
    return rename!(out, "flows" * "_$(string(f))" => :toflux)
end

function plotgroups!(df, ax, fun, group, xcol, ycol; kwargs...)
    for g in unique(df[!, group])
        foo = df[df[!, group] .== g, :]
        if isnothing(ycol)
            fun(ax, foo[!, xcol]; label = g, kwargs...)
        else
            fun(ax, foo[!, xcol], foo[!, ycol]; label = g, kwargs...)
        end
    end
end
