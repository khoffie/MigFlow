function estimate(model, model_args::NamedTuple, show_plt = true)
    ma = model_args
    mdl = model(ma)
    mles, preds = runoptim(mdl.mdl, mdl.lb, mdl.ub)
    out = format_mles(mles)

    net = calc_net_df(DataFrame(flows = ma.flows,
                                preds = preds,
                                fromdist = ma.fromdist,
                                todist = ma.todist))

    plt = [plotfit(ma.flows, preds),
           plotdist(ma.flows, preds, ma.dist),
           plotpop(ma.flows, preds, ma.frompop, ma.topop),
           plotnet(net)
           ## plotdens(ma[:flows], preds, ma[:fd], ma[:td])
           ]
    p = plot(plt..., plot_title = "LP $(round(out[end], digits = 0))",
             size = (1000, 600))
    if show_plt; display(p); end
    res = (out = out, net = net, preds = preds, mles = mles, plt = plt)
    return res
end

function instantiate(model, model_args)
    m = model(model_args)
    return m
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
            mles = @time(Turing.maximum_likelihood(mdl; lb = lb, ub = ub))
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
