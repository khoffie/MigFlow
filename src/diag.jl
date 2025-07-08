function calc_net(df, col)
    netf = combine(DataFrames.groupby(df, :fromdist), col => sum)
    rename!(netf, string(col) * "_sum" => :outflux)
    nett = combine(DataFrames.groupby(df, :todist), col => sum)
    rename!(nett, string(col) * "_sum" => :influx)
    net = innerjoin(netf, nett, on = [:fromdist => :todist])
    net.net = net.influx .- net.outflux
    net.total = net.influx .+ net.outflux
    net.nmr = net.net ./ net.total
    return net
end

function calc_net_df(df)
    net = calc_net(df, :flows)
    netp = calc_net(df, :preds)

    new = names(netp) .* "p"
    names(netp)
    rename!(netp, names(netp) .=> new)
    net = innerjoin(net, netp, on = :fromdist => :fromdistp)
    return net
end
