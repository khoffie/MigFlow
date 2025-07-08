function estimate(model; show_plt = true, optim_kwargs = (;))
    mles, preds = runoptim(mdl.mdl, mdl.lb, mdl.ub; optim_kwargs...)
    out = format_mles(mles)
    out = add_age_year(out, mdl.data)
    return diagchain(makechain(out), mdl)
end

function add_age_year(out, data)
    new_vals = vcat(out, [data.age, data.year])
    new_names = vcat(names(out)[1], ["age", "year"])
    return NamedArray(new_vals, (new_names,))
end

function format_mles(mles)
    out = NamedArray([mles.values.array..., mles.lp, Int(mles.optim_result.retcode)])
    nms = [string.(names(mles.values)...)..., "lp", "retcode"]
    setnames!(out, nms, 1)
    return out
end

function runoptim(mdl, lb, ub;
                  ad = ADTypes.AutoForwardDiff(),
                  inits = nothing,
                  reltol = nothing,
                  maxiters = nothing,
                  maxtime = nothing,
                  show_trace = false)
    attempt = 0
    while attempt < 5
        try
            mles = Turing.maximum_likelihood(mdl; lb = lb, ub = ub,
                                             adtype = ad,
                                             initial_params = inits,
                                             maxiters = maxiters,
                                             maxtime = maxtime,
                                             reltol = reltol,
                                             show_trace = show_trace,
                                             store_trace = false)
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

function makechain(out)
    coefs = (out[1 : end - 2])
    coefsnum = Float64.(coefs.array.array)
    nms = names(coefs)[1]
    return Chains(reshape(coefsnum, (1, :, 1)), nms)
end

## not needed because the lp Turing returns from optimization is
## actually the loglikelihood
# function makechain(out)
#     coefs = Float64.(collect(out.out[1 : end - 3]))
#     coefs = reshape(coefs, (1, length(coefs), 1))
#     return Chains(coefs, names(out.out)[1][1 : end - 3])
# end

# LL(mdl, out) = loglikelihood(mdl, makechain(out))
# LP(mdl, out) = logjoint(mdl, makechain(out))
