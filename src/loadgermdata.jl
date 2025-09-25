function load_data(a::String, y::Int, p::Float64, path::String;
                   only_positive::Bool, seed::Int = 123)
    di = CSV.read(joinpath(path, "districts.csv"), DataFrame)
    di = addlrd!(year(di, y))
    df = CSV.read(joinpath(path, "FlowDataGermans.csv"), DataFrame)
    df = year(age(df, a), y)
    Random.seed!(seed)
    if only_positive; df = pos(df); end
    if p < 1.0; df = sample_flows(df, p); end
    return (df = df, districts = di)
end

function sample_flows(flows::DataFrame, p::AbstractFloat)
    ## just to be save in case some years or age groups are passed
    ods = unique(flows, [:fromdist, :todist])[:, [:fromdist, :todist]]
    ## / 2 bc every pair is duplicated to make sure for o -> d the
    ## other direction is included as well
    nrows = Int(floor(nrow(ods) * p / 2))
    ods = ods[StatsBase.sample(1 : nrow(ods), nrows, replace = false), :]
    fromdist = vcat(ods.fromdist, ods.todist)
    todist = vcat(ods.todist, ods.fromdist)
    ## since for all o -> d, we also use all d -> o, and vice versa,
    ## we need to make sure we throw away all duplicated rows
    ods = unique(DataFrame(; fromdist, todist))
    return innerjoin(flows, ods, on = [:fromdist, :todist])
end

function addlrd!(districts::DataFrame)
    calc_lrd(x) = log.(x ./ median(x))
    districts.lrd = calc_lrd(districts.density)
    DataFrames.transform!(DataFrames.groupby(districts, :year),
                          :density => calc_lrd => :lrd)
    return districts
end

lrd(x) = log.(x ./ median(x))
radius(pop, dens) = sqrt.(pop ./ dens ./ 2π)

function sample_rows(df::DataFrame, p::AbstractFloat)
    nrows = nrow(df)
    n = Int(floor(nrows * p))
    df = df[Random.shuffle(1:nrows)[1 : n], : ]
    return(df)
end
