function chaindiag(chain, mdl, data)
    chn[:lp][end]
    p1 = plot(chn[:lp], label = "$(round(chn[:lp][end], digits = 2))")

    denscoefs = extract_coefs(chn[end, :, 1], "ζ")
    mat, p2 = plotdensrbf(denscoefs, mdl.data)

    geocoefs = extract_coefs(chn[end, :, 1], "η")
    dfgeo, p3 = plotgeorbf(geocoefs, mdl.data)

    df = DataFrame(;
                   data.df.fromdist,
                   data.df.todist,
                   data.df.flows,
                   preds = returned(mdl.mdl, chn[end])[1],
                   data.df.dist);

    p4 = plot(plotfit(df.flows, df.preds),
              title = "Cor: $(corround(log.(df.flows), log.(df.preds)))")

    p5 = plotdist(df, :preds)

    net = calc_net_df(df);
    p6 = plot(plotnet(net), title = "cor: $(corround(net.nmr, net.nmrp))")
    plts = [p1, p2, p3, p4, p5, p6]
    return (; df, net, plts)
end

function extract_coefs(chn, string)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end
