function gen_mdat(age::String, year::Int, p::AbstractFloat, path::String;
                  type::String, distscale::AbstractFloat = 100.0,
                  ndc::Signed, ngc::Signed, only_positive::Bool)

    data = load_data(age, year, p, path, only_positive);
    df = data.df
    districts = data.districts
    check_distcodes(df.fromdist, df.todist, districts.distcode)

    datafull = load_data(age, year, 1.0, path, only_positive);
    dffull = datafull.df

    if type == "joint"
        frompop = df.frompop
    elseif type == "conditional"
        out = combine(groupby(df, :fromdist), :flows => sum => :outflow)
        df2 = leftjoin(df, out, on = :fromdist)
        frompop = df2.outflow
    end

    dat = (
        flows = df.flows,
        ## outflow = df.outflow,
        fromdist = levelcode.(categorical(df.fromdist)),
        todist = levelcode.(categorical(df.todist)),
        fromdist_orig = df.fromdist,
        todist_orig = df.todist,
        agegroup = df.agegroup,
        year = df.year,
        frompop = frompop,
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
        ngc = ngc,
        radius = radius(districts.pop, districts.density),
        fpt = districts.pop,
        fromdistfull = levelcode.(categorical(dffull.fromdist)),
        topopfull = dffull.topop,
        distfull = dffull.dist
    )
    return dat
end

function check_distcodes(x, y, z)
    su(x) = sort(unique(x))
    @assert su(x) == su(y) "x and y are not equal"
    @assert su(y) == su(z) "y and z are not equal"
end
