using CSV, DataFrames, Turing, Mooncake, StatsBase, Random, Distributions,
    CategoricalArrays, NamedArrays, ADTypes, DynamicPPL, Serialization

include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/modelutils.jl")
include("../src/utils.jl")

include("../models/baseflow.jl")
include("../models/fundamental.jl")
include("../models/gravity.jl")
include("../models/norm.jl")

outp = "./output/"

function defmodel(m, age, year)
    data = load_data(age, year, 1.0, "../data/"; only_positive = true)
    return m == baseflow ? m(data; ndc = 16, ngcx = 5) : m(data)
end

function fitmodels(models, ages, years)
    for  m in models
        for a in ages
            name = "$(m)_$(a)"
            results = Vector{EstimationResult}(undef, length(years))
            Threads.@threads for i in eachindex(years)
                mdl = defmodel(m, a, years[i])
                results[i] = @time estimate(mdl)
            end
            serialize(joinpath(outp, name), results)
            println("$name done")
        end
    end
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
years = vcat(2000:2002, 2004:2017)
models =  [baseflow, fundamental, norm, gravity]

fitmodels([baseflow], ages, years)
