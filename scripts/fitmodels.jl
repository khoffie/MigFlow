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
    else
        mdl = m(data; trunc = trunc, norm = norm)
    end
    return mdl
end

function fitmodels(models, ages, years, trunc, norm, outp = "./output", prefit = false)
    if !isdir(outp); mkdir(outp); end
    for  m in models
        for a in ages

            name = modelname(m, trunc, norm, a)
            rs = Vector{EstimationResult}(undef, length(years))

            Threads.@threads for i in eachindex(years)
                if prefit
                    ## we prefit with 10% of rows, no truncation, no
                    ## normalization. Both cause problems with full
                    ## model
                    mdl = defmodel(m, a, years[i], .1, false, false)
                    out = estimate(mdl)
                    if norm
                        ## we need to smuggle in an init for beta
                        inits = vcat(out.ses.coef[1], 1.0, out.ses.coef[2:end])
                    else
                        inits = out.ses.coef
                    end
                    mdl = defmodel(m, a, years[i], 1.0, trunc, norm)
                    rs[i] = @time estimate(mdl, optim_kwargs = (; initial_params = inits, maxtime = 200))
                else
                    mdl = defmodel(m, a, years[i], 1.0, trunc, norm)
                    rs[i] = @time estimate(mdl, optim_kwargs = (; initial_params = inits, maxtime = 200))
                end
            end
            serialize(joinpath(outp, name), rs)
            println("$name done")
        end
    end
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
years = vcat(2000:2002, 2004:2017)
## models =  [baseflow, fundamental, norm, gravity, baseflownormalized]

# fitmodels([fundamental, gravity], ages, years, true, false, "./output", false)
# fitmodels([fundamental], ages, years, true, true, "./output", false)
## fitmodels([baseflow], ages, years, true, false, "./output", true)
fitmodels([baseflow], ages, years, true, true, "./output", true)
