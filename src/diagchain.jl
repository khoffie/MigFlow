function chaindiag(chn, mdl, data)
    p1 = plot(chn[:lp])

    m = argmax(chn[:lp].data)
    chn = chn[m[1], :, m[2]]
    df = DataFrame(;
                   data.df.fromdist,
                   data.df.todist,
                   data.df.flows,
                   preds = returned(mdl.mdl, chn)[1],
                   data.df.dist);

    p2 = plot(plotfit(df.flows, df.preds),
              title = "Cor: $(corround(log.(df.flows), log.(df.preds)))")

    if occursin("dist", names(df)[1])
        p3 = plotdist(df, :preds)
    else
        p3 = nothing
    end

    net = calc_net_df(df);
    p4 = plot(plotnet(net), title = "cor: $(corround(net.nmr, net.nmrp))")

    denscoefs = extract_coefs(chn[end, :, 1], "ζ")
    if length(denscoefs) > 0
        mat, p5 = plotdensrbf(denscoefs, mdl.data)
    else
        mat, p5 = nothing, nothing
    end

    geocoefs = extract_coefs(chn[end, :, 1], "η")
    if length(geocoefs) > 0
        geo, p6 = plotgeorbf(geocoefs, mdl.data)
    else
        geo, p6 = nothing, nothing
    end

    plts = [p1, p2, p3, p4, p5, p6]
    return (; df, net, mat, geo, plts)
end

function extract_coefs(chn, string)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end
