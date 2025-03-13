function read_flows(file::String, p::AbstractFloat,
                    ages, years, pos_only = true)
    df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
    df = sample_flows(df, p)
    check_age(age) = age in ages
    check_year(year) = year in years
    df = filter(:agegroup => check_age, df)
    df = filter(:year => check_year, df)
    if pos_only; df = df[df.flows .> 0, :]; end
    return df
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

function add_lrd(districts::DataFrame)
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

# function load_data(datapath::String,
#                    p::AbstractFloat = 1.0,
#                    age::String = "30-50",
#                    year::Signed = 2017,
#                    seed::Int = 1234)
#     Random.seed!(seed)
#     germ = loadallGermData(datapath; p = p, positive_only = true)
#     fl = germ.flows[germ.flows.year .== year, :]
#     fl = fl[fl.agegroup .== age, :]
#     return fl, germ.geog
# end

# function loadallGermData(path = "data"; p::AbstractFloat,
#                          positive_only::Bool = true)
#     geog = loadGermGeog(path)
#     geog = sample_germ(geog, p)
#     flows = loadGermFlows(path, levels(geog.distcode), positive_only)
#     return (geog=geog, flows=flows)
# end

# function loadGermGeog(path)
#     geog = CSV.read(joinpath(path, "districts.csv"), DataFrame)
#     geog.distcode = categorical(geog.distcode)
#     sort!(geog,[:distcode])
#     meddens = median(geog.density)
#     geog.logreldens = log.(geog.density ./ meddens)

#     geog
# end

# function sample_germ(geog::DataFrame, p::AbstractFloat)
#     geog = sample_rows(geog, p)
#     CategoricalArrays.droplevels!(geog.distcode)
#     return geog
# end

# function loadGermFlows(path, distcodes, positive_only)
#     fl = CSV.read(joinpath(path, "FlowDataGermans.csv"), DataFrame)
#     fl = filter(row -> row.fromdist in distcodes, fl)
#     fl = filter(row -> row.todist in distcodes, fl)
#     fl.fromdist = categorical(fl.fromdist, levels = distcodes)
#     fl.todist = categorical(fl.todist, levels = distcodes)
#     if positive_only
#         fl = filter(row -> row.flows .> 0, fl)
#     end
#     return fl
# end
