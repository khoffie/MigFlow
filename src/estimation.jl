function estimate(mdl; show_plt = true, optim_kwargs = (;))
    mles = runoptim(mdl.mdl, mdl.lb, mdl.ub; optim_kwargs...)
    out = format_mles(mles)
    out = add_age_year(out, mdl.data)
    return (chn = makechain(out), mdl = mdl)
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
            mles = Turing.maximum_a_posteriori(mdl; lb = lb, ub = ub,
                                             adtype = ad,
                                             initial_params = inits,
                                             maxiters = maxiters,
                                             maxtime = maxtime,
                                             reltol = reltol,
                                             show_trace = show_trace,
                                             store_trace = false)
##            mles = @time(Turing.maximum_a_posteriori(mdl; lb = lb, ub = ub))
            return mles
        catch e
            if e isa DomainError
                println("Domainerror, retrying... Attempt $attempt")
                attempt += 1
            else
                rethrow(e)
            end
        end
    end
    error("Optimization failed after $attempt attempts")
end

function makechain(out)
    coefs = vcat(out[1 : end - 2], recodeage(String(out["age"])), out["year"])
    coefsnum = Float64.(coefs.array.array)
    nms = names(out)[1]
    return Chains(reshape(coefsnum, (1, :, 1)), nms)
end

function recodeage(a::String)
    a == "below18" && return 1
    a == "18-25" && return 2
    a == "25-30" && return 3
    a == "30-50" && return 4
    a == "50-65" && return 5
    a == "above65" && return 6
end

function recodeage(i::Int)
    i == 1 && return "below18"
    i == 2 && return "18-25"
    i == 3 && return "25-30"
    i == 4 && return "30-50"
    i == 5 && return "50-65"
    i == 6 && return "above65"
end
