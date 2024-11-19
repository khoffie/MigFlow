sample_germ = function(geog; p = .1)
    geog = sample_rows(geog, p)
    CategoricalArrays.droplevels!(geog.distcode)
    return geog
end

sample_us = function(geog, flows; p = .1)
    geog = sample_rows(geog, p)
    CategoricalArrays.droplevels!(geog.countyid)
    countycodes = unique(geog.countyid)
    flows = filter(row -> row.fromcounty in countycodes, flows)
    flows = filter(row -> row.tocounty in countycodes, flows)
    return (geog, flows)
end

sample_rows = function(df, p)
    nrows = nrow(df)
    n = Int(floor(nrows * p))
    df = df[Random.shuffle(1:nrows)[1 : n], : ]
    return(df)
end

