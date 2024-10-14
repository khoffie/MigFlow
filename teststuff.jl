using StatsBase
germ = loadallGermData(; sample = false)
us = loadallUSdata()

sample_rows = function(df, p)
    nrows = nrow(df)
    n = Int(floor(nrows * p))
    df = df[Random.shuffle(1:nrows)[1 : n], : ]
    return(df)
end

sample_rows_germ = function(df; p = .1)
    df = sample_rows(df, p = p)
    CategoricalArrays.droplevels!(df.distcode)
    return df
end

sample_rows_us = function(df; p = .1)
    df = sample_rows(df, p = p)
    return df
end

    

test_germ = sample_rows_germ(germ.geog)
test_us = sample_rows_us(us.geog)

