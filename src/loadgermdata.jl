function load_data(a::String, y::Int, p::Float64, path::String;
                   positive::Bool, full::Bool)
    di = CSV.read(joinpath(path, "districts.csv"), DataFrame)
    di = addlrd!(di)
    dffull = CSV.read(joinpath(path, "FlowDataGermans.csv"), DataFrame)
    dffull = year(age(dffull, a), y)
    if positive; dffull = pos(dffull); end
    if p < 1.0
        df = sample_flows(dffull, p)
    elseif p == 1.0
        df = dffull
    end
    if !full
        out = (df = df, districts = year(di, y))
    elseif full
        dffull = dffull[!, [:fromdist, :todist, :dist]]
        out = (df = df, districts = year(di, y), dffull)
    end
    return out
end

function sample_flows(flows::DataFrame, p::AbstractFloat)
    ods = unique(flows, [:fromdist, :todist])[:, [:fromdist, :todist]]
    ## / 2 bc every pair is duplicated to make sure for o -> d the
    ## other direction is included as well
    nrows = Int(floor(nrow(ods) * p / 2))
    ods = ods[StatsBase.sample(1 : nrow(ods), nrows, replace = false), :]
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

function joinlrd(df, districts)
    di = select(districts, :distcode, :year, :lrd => :fromdens)
    leftjoin!(df, di, on = [:fromdist => :distcode, :year])
    di = select(districts, :distcode, :year, :lrd => :todens)
    leftjoin!(df, di, on = [:todist => :distcode, :year])
    return df
end
