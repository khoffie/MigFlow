function netdf(r::EstimationResult, percent = true)
    return addmeta(r, addnetcols(calcnet(modeldf(r)), percent))
end

function netdf(df::DataFrame, percent = true)
    return addnetcols(calcnet(df, Flows()), percent)
end

function calcnet(df::DataFrame)
    dfout = combine(DataFrames.groupby(df, [:fromdist]),
                    :flows => sum => :outflux,
                    :preds => sum => :outfluxp)
    dfin = combine(DataFrames.groupby(df, [:todist]),
                   :flows => sum => :influx,
                   :preds => sum => :influxp)
    net = innerjoin(dfout, dfin, on = :fromdist => :todist)
    net = innerjoin(net, unique(df, :fromdist)[!, [:fromdist, :A]], on = :fromdist)
    return rename!(net, :fromdist => :lc)
end

function calcnet(df::DataFrame, type::Flows)
    df2 = calcnet(df)
    df2.agegroup .= unique(df.agegroup)[1]
    df2.year .= unique(df.year)[1]
    return df2
end

function addnetcols(df::DataFrame)
    net = df
    net.net = net.influx .- net.outflux
    net.total = net.influx .+ net.outflux
    net.asyma = net.net ./ net.total
    net.nmra = net.net ./ net.A
    net.netp = net.influxp .- net.outfluxp
    net.totalp = net.influxp .+ net.outfluxp
    net.asymap = net.netp ./ net.totalp
    net.nmrap = net.netp ./ net.A
    net.diff = net.nmra .- net.nmrap
    return net
end
