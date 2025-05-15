# function read_flows(file::String, p::AbstractFloat,
#                     ages, years, pos_only = true)
#     df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
#     df = sample_flows(df, p)
#     check_age(age) = age in ages
#     check_year(year) = year in years
#     df = filter(:agegroup => check_age, df)
#     df = filter(:year => check_year, df)
#     if pos_only; df = df[df.flows .> 0, :]; end
#     return df
# end

function load_data(a::String, y::Int, p::Float64, path::String)
    di = CSV.read(joinpath(path, "districts.csv"), DataFrame)
    di = addlrd!(di)
    df = CSV.read(joinpath(path, "FlowDataGermans.csv"), DataFrame)
    df = year(age(pos(df), a), y)
    df = sample_flows(df, p)
    df = joinlrd(df, di)
    return (df = df, districts = year(di, y))
end

function sample_flows(flows::DataFrame, p::AbstractFloat)
    ods = unique(flows, [:fromdist, :todist])[:, [:fromdist, :todist]]
    ## / 2 bc every pair is duplicated to make sure for o -> d the
    ## other direction is included as well
    nrows = Int(floor(nrow(ods) * p / 2))
    ods = ods[StatsBase.sample(1 : nrow(ods), nrows), :]
    fromdist = vcat(ods.fromdist, ods.todist)
    todist = vcat(ods.todist, ods.fromdist)
    ods = DataFrame(; fromdist, todist)
    return innerjoin(flows, ods, on = [:fromdist, :todist])
end

function addlrd!(districts::DataFrame)
    calc_lrd(x) = log.(x ./ median(x))
    districts.lrd = calc_lrd(districts.density)
    DataFrames.transform!(groupby(districts, :year),
                          :density => calc_lrd => :lrd)
    return districts
end

function joinlrd(flows, districts)
    flows = innerjoin(flows, districts[:, [:distcode, :lrd, :year]],
                      on = [:fromdist => :distcode, :year])
    rename!(flows, :lrd => :fromdens)
    flows = innerjoin(flows, districts[:, [:distcode, :lrd, :year]],
                      on = [:todist => :distcode, :year])
    rename!(flows, :lrd => :todens)
    return flows
end
