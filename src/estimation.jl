function estimate(model, model_args::NamedTuple, show_plt = true)
    ma = model_args
    mdl = model(ma)
    mles, preds = runoptim(mdl.mdl, mdl.lb, mdl.ub)
    out = format_mles(mles)

    net = calc_net_df(DataFrame(flows = ma.flows,
                                preds = preds,
                                fromdist = ma.fromdist,
                                todist = ma.todist))

    densdesir, pdens = evaldens(out, ma)
    geo, pgeo = evalgeo(out, ma)
    df = DataFrame(flows = ma.flows, preds = preds, dist = ma.dist)
    plt = [
        plotfit(ma.flows, preds),
        plotdist(df, :preds, 100),
        plotpop(ma.flows, preds, ma.frompop, ma.topop),
        plotnet(net),
        pdens,
        pgeo
#           plotdens(ma[:flows], preds, ma[:fd], ma[:td])
           ]
    p = plot(plt[1 : 4]..., plot_title = "LP $(round(out[end], digits = 0))",
             size = (1000, 600))
    if show_plt; display(p); end
    res = (out = out, net = net, preds = preds,
           dens = densdesir, geo = geo, mles = mles, plt = plt)
    return res
end

function evaldens(out, d)
    if any(occursin.("kd", names(out)[1]))
        kds = out[["kd[$i]" for i in 1 : d.ndc]]
        return evaldensitycheby(kds, d.dmin, d.dmax)
    else
        return nothing, nothing
    end
end

function evalgeo(out, d)
    if any(occursin.("kg", names(out)[1]))
        kgs = out[["kg[$i]" for i in 1 : d.ngc]]
        return evalgeocheby(kgs, d.distcode, d.xcoord, d.ycoord, true)
    else
        return nothing, nothing
    end
end

function format_mles(mles)
    out = NamedArray([mles.values.array..., mles.lp])
    nms = [string.(names(mles.values)...)..., "lp"]
    setnames!(out, nms, 1)
    return out
end

function runoptim(mdl, lb, ub)
    attempt = 0
    while attempt < 5
        try
            mles = Turing.maximum_likelihood(mdl; lb = lb, ub = ub)
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
