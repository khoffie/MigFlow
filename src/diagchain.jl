function diagchain(chn, mdl)
    p1 = plot(chn[:lp])

    m = argmax(chn[:lp].data)
    chn = chn[m[1], :, m[2]]
    df = DataFrame(
                   fromdist = mdl.mdl.args.from,
                   todist = mdl.mdl.args.to,
                   flows = mdl.mdl.args.Y,
                   preds = returned(mdl.mdl, chn)[1],
                   dist = mdl.mdl.args.D);
    dev = round(deviance(df.flows, df.preds), digits = 2)
    p2 = plot(plotfit(df.flows, df.preds),
              title = "Mean deviance: $(dev)")

    if occursin("dist", names(df)[1])
        p3 = plotdist(df, :preds)
    else
        p3 = nothing
    end

    net = calc_net_df(df);
    p4 = plot(plotnet(net), title = "cor: $(corround(net.nmr, net.nmrp))")

    denscoefs = extract_coefs(chn[end, :, 1], "ζ")
    if length(denscoefs) > 0
        mat, p5 = plotdensrbf(denscoefs, mdl.mdl.args)
    else
        mat, p5 = nothing, nothing
    end

    geocoefs = extract_coefs(chn[end, :, 1], "η")
    if length(geocoefs) > 0
        geo, p6 = plotgeorbf(geocoefs, mdl.mdl.args)
    else
        geo, p6 = nothing, nothing
    end

    plts = [p1, p2, p3, p4, p5, p6]
    return (; df, net, mat, geo, plts, chn)
end

function extract_coefs(chn, string)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end

function extract_sample(chn, type = "best")
    if type == "best"
        m = argmax(chn[:lp].data)
        chn = chn[m[1], :, m[2]]
    end
    return chn
end

function deviance(y, p)
    loss = zeros(length(y))
    for i in eachindex(y)
        loss[i] = y[i] * log(y[i] / p[i]) - (y[i] - p[i])
    end
    return 2mean(loss)
end

corround(x, y) = round(cor(x, y), digits = 2)
