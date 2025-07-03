function chaindiag(chain, mdl, data)
    df = DataFrame(;
                   data.df.fromdist,
                   data.df.todist,
                   data.df.flows,
                   preds = returned(mdl.mdl, chn[end])[1],
                   data.df.dist);

    chn[:lp][end]
    p1 = plot(chn[:lp], label = "$(round(chn[:lp][end], digits = 2))")
    p2 = plot(plotfit(df.flows, df.preds),
              title = "Cor: $(corround(log.(df.flows), log.(df.preds)))")
    p3 = plotdist(df, :preds)
    net = calc_net_df(df);
    p4 = plot(plotnet(net), title = "cor: $(corround(net.nmr, net.nmrp))")
    denscoefs = extract_coefs(chn[end, :, 1], "ζ")
    mat, p5 = plotdensrbf(denscoefs, mdl.data)
    geocoefs = extract_coefs(chn[end, :, 1], "η")
    dfgeo, p6 = plotgeorbf(geocoefs, mdl.data)

    plts = [p1, p2, p3, p4, p5, p6]
    return (; df, net, mat, plts)
end

function extract_coefs(chn, string)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end
