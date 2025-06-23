function estimate(model, data::NamedTuple;
                  ad = ADTypes.AutoForwardDiff(), show_plt = true,
                  kwargs...)
    mdl = model(data; kwargs...)
    mles, preds = runoptim(mdl.mdl, mdl.lb, mdl.ub, ad)
    out = format_mles(mles)
    data = mdl.data

    out = NamedArray(out ∪ [data.age, data.year],
                     names(out)[1] ∪ ["age", "year"])

    net = calc_net_df(DataFrame(flows = data.Y,
                                preds = preds,
                                fromdist = data.from,
                                todist = data.to))

    densdesir, pdens = evaldens(out, data)
    geo, pgeo = evalgeo(out, data)
    df = DataFrame(flows = data.Y, preds = preds, dist = data.D)
    plt = [
        plotfit(data.Y, preds),
        plotdist(df, :preds),
        plotpop(data.Y, preds, data.A[data.from], data.P[data.to]),
        plotnet(net),
        pdens,
        pgeo
#           plotdens(ma[:flows], preds, ma[:fd], ma[:td])
           ]
    p = plot(plt[1 : 4]..., plot_title = "LP $(round(out["lp"], digits = 0))",
             size = (1000, 600))
    if show_plt; display(p); end
    res = (out = out, net = net, preds = preds,
           dens = densdesir, geo = geo, mles = mles, plt = plt)
    display(res.out)
    println("")
    return res
end

function format_mles(mles)
    out = NamedArray([mles.values.array..., mles.lp])
    nms = [string.(names(mles.values)...)..., "lp"]
    setnames!(out, nms, 1)
    return out
end

function runoptim(mdl, lb, ub, ad)
    attempt = 0
    while attempt < 5
        try
            mles = Turing.maximum_likelihood(mdl; lb = lb, ub = ub,
                                             adtype = ad)
##            mles = @time(Turing.maximum_a_posteriori(mdl; lb = lb, ub = ub))
            return mles, predict(mdl, mles)
        catch e
            if e isa DomainError
                println("Domainerror, retrying... Attempt $attempt")
                attempt += 1
            else
                rethrow(e)
            end
        end
    end
    error("Optimization failed after $attempt attempts")  # Ensure a clear failure message
end

function predict(mdl, mles)
    ## mdl is an instantiated Turing model
    params = string.(names(mles.values)...)
    preds = generated_quantities(mdl, mles.values, params)
    return preds
end
