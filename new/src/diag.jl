function calc_net(flows, col)
    netf = combine(groupby(flows, :fromdist), col => sum)
    rename!(netf, string(col) * "_sum" => :influx)
    nett = combine(groupby(flows, :todist), col => sum)
    rename!(nett, string(col) * "_sum" => :outflux)
    net = innerjoin(netf, nett, on = [:fromdist => :todist])
    net.net = net.influx .- net.outflux
    net.total = net.influx .+ net.outflux
    net.nmr = net.net ./ net.total
    return net
end

function calc_net_df(flows)
    net = calc_net(flows, :flows)
    netp = calc_net(flows, :preds)

    new = names(netp) .* "p"
    names(netp)
    rename!(netp, names(netp) .=> new)
    net = innerjoin(net, netp, on = :fromdist => :fromdistp)
    return net
end

function evalgeocheby(coefs, districts,
                      show_plt = false)
    xmin, xmax = minmax(districts.xcoord)
    ymin, ymax = minmax(districts.ycoord)
    geocheby = defgeocheby(coefs, xmin, xmax, ymin, ymax)
    geo = geocheby.(districts.xcoord, districts.ycoord)
    ##geo = ForwardDiff.value.(geo)
    df = DataFrame(; districts.distcode,
                   districts.xcoord,
                   districts.ycoord, geo)

    ratio = (ymax - ymin) / (xmax - xmin)
    width = 600
    p = scatter(df.xcoord, df.ycoord,
                marker_z = df.geo,
                markersize = 10, size = (600, width * ratio))
    if show_plt; display(p); end
    return df, p
end

function evaldensitycheby(coefs, densmin, densmax,
                          vals = range(densmin, densmax, 100))
    densitycheby = defdensitycheby(coefs, densmin, densmax)
    df = ((fromdens = fd,
          todens = td,
          funval = densitycheby(fd, td)) for fd in vals, td in vals)
    df = DataFrame(df)
    p = heatmap(vals, vals, reshape(df.funval, (100, 100)))
    display(p)
    return df, p
end
