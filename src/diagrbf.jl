function plotgeorbf(n::NamedTuple)
    geocoefs = extract_coefs(n.chn[end, :, 1], "η")
    data = n.mdl.mdl.args
    y = n.chn[:year][1]
    return plotgeorbf(geocoefs, data)
end

function plotgeorbf(coefs::Vector{Float64},
                    data::NamedTuple, title)
    xcoord = data.xcoord
    ycoord = data.ycoord
    cx = data.cxgeo
    cy = data.cygeo
    scale = data.geo_scale
    return plotgeorbf_(coefs, xcoord, ycoord, cx, cy, scale)
end

function plotgeorbf_(coefs::Vector{Float64},
                     xcoord::Vector{Float64},
                     ycoord::Vector{Float64},
                     cx::Vector{Float64},
                     cy::Vector{Float64},
                     scale::Float64,
                     title::String)
    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    geo = [interpolant(rbf, xcoord[i], ycoord[i],
                       coefmat(coefs ./ 10, length(cx), length(cy)),
                       cx, cy, scale) for i in eachindex(xcoord)];
    geo = geo .- mean(geo)
    dfgeo = DataFrame(; xcoord, ycoord, geo)
    ratio = (ymax - ymin) / (xmax - xmin)
    ratio = 1.33
    width = 600
    p = scatter(xcoord, ycoord,
                marker_z = dfgeo.geo,
                markersize = 10, size = (600, width * ratio),
                label = "", yflip = false, title = title)
    res = collect(IterTools.product(cx, cy))
    xvals = getindex.(res, 1)
    yvals = getindex.(res, 2)
    scatter!(vec(xvals), vec(yvals), color = :red)
    return dfgeo, p
end

function plotdensrbf(coefs::Vector{Float64}, data::NamedTuple)
    cx = data.cx
    cy = data.cy
    scale = data.rbf_scale
    R = data.R
    return plotdensrbf_(coefs, cx, cy, scale, R)
end

function plotdensrbf_(coefs::Vector{Float64},
                      cx::Vector{Float64},
                      cy::Vector{Float64},
                      scale::Float64,
                      R::Vector{Float64})
    Rmin, Rmax = extrema(R)
    vals = range(Rmin, Rmax, 1000)
    s = Int(sqrt(length(coefs)))
    mat = [interpolant(rbf, Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale)
           for Rto in vals, Rfrom in vals];
    mat = mat .- mean(mat)
    p = heatmap(mat)
    return mat, p
end
