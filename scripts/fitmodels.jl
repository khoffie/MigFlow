using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, ADTypes, DynamicPPL, Serialization,
    SpecialFunctions, LogExpFunctions, Statistics

include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/modelutils.jl")
include("../src/utils.jl")

include("../models/TruncatedPoisson.jl")
include("../models/baseflow.jl")
include("../models/fundamental.jl")
include("../models/gravity.jl")

outp = "./output/"

function defmodel(m, age, year, p, trunc, norm)
    data = load_data(age, year, p, "../data/"; only_positive = true)
    if m == baseflow
        mdl = m(data; ndc = 16, ngcx = 5, trunc = trunc, norm = norm)
    elseif m == fundamental
        mdl = m(data; trunc = trunc, norm = norm)
    end
    return mdl
end

function fitmodels(models, ages, years, trunc, norm, outp = "./output",
                   prefit = false, suffix = nothing)
    if !isdir(outp); mkdir(outp); end
    for  m in models
        for a in ages

            name = "$(m)_$(a)"
            if !trunc; name * "_nontruncated"; end
            if norm; name * "normalized"; end
            if !isnothing(suffix); name = name * "_suffix"; end

            results = Vector{EstimationResult}(undef, length(years))
            Threads.@threads for i in eachindex(years)
                if prefit
                    mdl = defmodel(m, a, years[i], .1, false, false)
                    out = estimate(mdl)
                    ## we need to smuggle in an init for beta
                    inits = vcat(out.ses.coef[1], 1.0, out.ses.coef[2:end])
                    mdl = defmodel(m, a, years[i], 1.0, trunc, norm)
                    results[i] = @time estimate(mdl, optim_kwargs = (; inits = inits))
                else
                    mdl = defmodel(m, a, years[i], 1.0, trunc, norm)
                    results[i] = @time estimate(mdl)
                end
            end
            serialize(joinpath(outp, name), results)
            println("$name done")
        end
    end
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
years = vcat(2000:2002, 2004:2017)
## models =  [baseflow, fundamental, norm, gravity, baseflownormalized]

fitmodels([fundamental], ages, years, true, false, "./output", false)
fitmodels([fundamental], ages, years, true, true, "./output", false)
