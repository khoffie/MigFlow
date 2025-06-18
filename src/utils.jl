logistic(x) = 1 / (1 + exp(-x))

pos(df) = df[df.flows .> 0.0, :]
age(df, age) = filter(:agegroup => n -> n == age, df)
year(df, y) = filter(:year => n -> n == y, df)
lrd(x) = log.(x ./ median(x))
radius(pop, dens) = sqrt.(pop ./ dens ./ 2π)

function sample_rows(df::DataFrame, p::AbstractFloat)
    nrows = nrow(df)
    n = Int(floor(nrows * p))
    df = df[Random.shuffle(1:nrows)[1 : n], : ]
    return(df)
end
