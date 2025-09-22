struct ModelWrapper
    mdl::Model # Turing Model
    lb::Vector{Float64} # Lower bound for optimization
    ub::Vector{Float64} # Upper bound for optimization
    data::NamedTuple # Meta data of fit: Age and Year
end

struct EstimationResult
    chn::Chains # Turing MCMC Chain, optimized parameters are stored here
    mdl::ModelWrapper
    prd::Vector{Float64} # Model predictions
    maps::Turing.Optimisation.ModeResult
end

function estimate(mdl::ModelWrapper; show_plt = true, optim_kwargs = (;))
    maps = runoptim(mdl; optim_kwargs...)
    chn = makechain(maps)
    prd = Turing.returned(mdl.mdl, chn)[1]
    return EstimationResult(chn, mdl, prd, maps)
end

function runoptim(mdl::ModelWrapper;
                  ad = ADTypes.AutoMooncake(),
                  inits = nothing,
                  reltol = nothing,
                  maxiters = nothing,
                  maxtime = nothing,
                  show_trace = false)
    attempt = 0
    while attempt < 5
        try
            mles = Turing.maximum_a_posteriori(mdl.mdl; lb = mdl.lb, ub = mdl.ub,
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

function makechain(maps::Turing.Optimisation.ModeResult)
    estims = maps.values.array
    nms  = names(maps.values)[1]
    return Chains(reshape(estims, (1, :, 1)), nms)
end
