logistic(x) = 1 / (1 + exp(-x))

function sample_rows(df::DataFrame, p::AbstractFloat)
    nrows = nrow(df)
    n = Int(floor(nrows * p))
    df = df[Random.shuffle(1:nrows)[1 : n], : ]
    return(df)
end

function gen_mdat(df::DataFrame, districts::DataFrame;
                  distscale::AbstractFloat = 100.0,
                  ndc::Signed, ngc::Signed)
    districts = unique(districts, :distcode)
    check_distcodes(df.fromdist, df.todist, districts.distcode)
    data = (
        flows = df.flows,
        fromdist = levelcode.(categorical(df.fromdist)),
        todist = levelcode.(categorical(df.todist)),
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
    return data
end

function check_distcodes(x, y, z)
    su(x) = sort(unique(x))
    @assert su(x) == su(y) "x and y are not equal"
    @assert su(y) == su(z) "y and z are not equal"
end
