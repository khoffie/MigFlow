 Base.@kwdef mutable struct MetaData
     ## mutable because age and year are assigned before fitting and
     ## lp can only be assigned afterwards
     model:: String # model, on of [baseflow, fundamental, normalized, gravity]
     age::String
     year::Int
     lp::Union{Nothing,Float64} = nothing
end

struct ModelWrapper
    mdl::DynamicPPL.Model # Turing Model
    lb::Vector{Float64} # Lower bound for optimization
    ub::Vector{Float64} # Upper bound for optimization
    meta::MetaData # Meta data of fit: Age, Year, LP
end

struct EstimationResult
    chn::Chains # Turing MCMC Chain, optimized parameters are stored
                # here, because then postprocessing functions like
                # plotting work with optimization and sampling
    mdl::ModelWrapper
    prd::Vector{Float64} # Model predictions
    ses::DataFrame ## standard errors
    cvm::NamedMatrix ## covariance matrix
end

function estimate(mdl::ModelWrapper, method=Optim.LBFGS(); ret_maps = false, optim_kwargs = (;))
    maps = runoptim(mdl, method; optim_kwargs...)
    mdl.meta.lp = maps.lp
    chn = makechain(maps)
    prd = Turing.returned(mdl.mdl, chn)[1]
    ses = setable(maps)
    cvm = vcov(maps)
    res = EstimationResult(chn, mdl, prd, ses, cvm)
    if ret_maps; res = res, maps; end
    return res
end

function runoptim(mdl::ModelWrapper, method;
                  ad = ADTypes.AutoMooncake(),
                  inits = nothing,
                  reltol = nothing, abstol = nothing,
                  maxiters = nothing,
                  maxtime = nothing,
                  show_trace = false)
    attempt = 0
    while attempt < 5
        try
            mles = Turing.maximum_a_posteriori(mdl.mdl, method; lb = mdl.lb, ub = mdl.ub,
                                             adtype = ad,
                                             initial_params = inits,
                                             maxiters = maxiters,
                                             maxtime = maxtime,
                                             reltol = reltol, abstol=abstol,
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

function setable(maps)
    ses = DataFrame(coeftable(maps))[!, 1:3]
    ses = rename(ses, ["Name", "Coef.", "Std. Error"] .=> ["name", "coef", "se"])
    ses.name = replace.(ses.name, "_raw" => "", "[" => "", "]" => "")
    return ses
end
