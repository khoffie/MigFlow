mm(v) = minimum(v), maximum(v)
logistic(x) = 1 / (1 + exp(-x))

pos(df) = df[df.flows .> 0.0, :]
age(df, age) = filter(:agegroup => n -> n == age, df)
year(df, y) = filter(:year => n -> n == y, df)

lrd(x) = log.(x ./ median(x))
addlrd(df) = (df.lrd = combine(groupby(df, :year), :density => lrd)[!, 2])

function joinlrd(df, districts)
    di = select(districts, :distcode, :year, :lrd => :fromdens)
    leftjoin!(df, di, on = [:fromdist => :distcode, :year])
    di = select(districts, :distcode, :year, :lrd => :todens)
    leftjoin!(df, di, on = [:todist => :distcode, :year])
    return df
end

function sample_rows(df::DataFrame, p::AbstractFloat)
    nrows = nrow(df)
    n = Int(floor(nrows * p))
    df = df[Random.shuffle(1:nrows)[1 : n], : ]
    return(df)
end

function gen_mdat(data::NamedTuple;
                  distscale::AbstractFloat = 100.0,
                  ndc::Signed, ngc::Signed)
    df = data.df
    districts = data.districts
    districts = unique(districts, :distcode)
    check_distcodes(df.fromdist, df.todist, districts.distcode)
    dat = (
        flows = df.flows,
        fromdist = levelcode.(categorical(df.fromdist)),
        todist = levelcode.(categorical(df.todist)),
        fromdist_orig = df.fromdist,
        todist_orig = df.todist,
        agegroup = df.agegroup,
        year = df.year,
        frompop = df.frompop,
        topop = df.topop,
        dist = df.dist,
        fromdens = df.fromdens,
        todens = df.todens,
        distcode = districts.distcode,
        xcoord = districts.xcoord,
        ycoord = districts.ycoord,
        dmin = minimum(districts.lrd),
        dmax = maximum(districts.lrd),
        distscale = distscale,
        ndc = ndc,
        ngc = ngc
    )
    return dat
end

function check_distcodes(x, y, z)
    su(x) = sort(unique(x))
    @assert su(x) == su(y) "x and y are not equal"
    @assert su(y) == su(z) "y and z are not equal"
end
