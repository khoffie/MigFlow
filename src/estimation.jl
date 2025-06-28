function estimate(model, data::NamedTuple;
                  show_plt = true,
                  model_kwargs = (;),
                  optim_kwargs = (;))
    mdl = model(data; model_kwargs...)
    mles, preds = runoptim(mdl.mdl, mdl.lb, mdl.ub; optim_kwargs...)
    out = format_mles(mles)
    data = mdl.data
    out = add_age_year(out, data)

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
           ]
    p = plot(plt[1 : 4]..., plot_title = "LP $(round(out["lp"], digits = 0))",
             size = (1000, 600))
    if show_plt; display(p); end
    res = (out = out, net = net, preds = preds,
           dens = densdesir, geo = geo, mles = mles, plt = plt)
##    display(res.out)
    println("")
    return res
end

function add_age_year(out, data)
    new_vals = vcat(out, [data.age, data.year])
    new_names = vcat(names(out)[1], ["age", "year"])
    return NamedArray(new_vals, (new_names,))
end

function format_mles(mles)
    out = NamedArray([mles.values.array..., mles.lp])
    nms = [string.(names(mles.values)...)..., "lp"]
    setnames!(out, nms, 1)
    return out
end

function runoptim(mdl, lb, ub;
                  ad = ADTypes.AutoForwardDiff(),
                  inits = nothing,
                  reltol = nothing,
                  maxiters = nothing,
                  show_trace = false)
    attempt = 0
    while attempt < 5
        try
            mles = Turing.maximum_likelihood(mdl; lb = lb, ub = ub,
                                             adtype = ad,
                                             initial_params = inits,
                                             maxiters = maxiters,
                                             reltol = reltol,
                                             show_trace = show_trace)
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
    vs = Tuple(mles.values)
    ks = Tuple(names(mles.values)[1])
    return returned(mdl, NamedTuple{ks}(vs))
end
