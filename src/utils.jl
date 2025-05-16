mm(v) = minimum(v), maximum(v)
logistic(x) = 1 / (1 + exp(-x))

pos(df) = df[df.flows .> 0.0, :]
age(df, age) = filter(:agegroup => n -> n == age, df)
year(df, y) = filter(:year => n -> n == y, df)
lrd(x) = log.(x ./ median(x))
## Should be save because even if population data is bad, in just
## about any data source density will be calculated as pop / area
radius(pop, dens) = sqrt.(pop ./ dens ./ 2π)

function sample_rows(df::DataFrame, p::AbstractFloat)
    nrows = nrow(df)
    n = Int(floor(nrows * p))
    df = df[Random.shuffle(1:nrows)[1 : n], : ]
    return(df)
end
