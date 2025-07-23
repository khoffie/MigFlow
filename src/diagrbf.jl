function plotgeorbf(n::NamedTuple, shp::DataFrame, clim = nothing)
    geocoefs = extract_coefs(n.chn[end, :, 1], "η")
    data = n.mdl.mdl.args
    y = n.chn[:year][1]
    return plotgeorbf(geocoefs, data, shp, "Year $y", clim)
end

function plotgeorbf(coefs::Vector{Float64},
                    data::NamedTuple, shp::DataFrame, title, clim)
    xcoord = data.xcoord
    ycoord = data.ycoord
    cx = data.cxgeo
    cy = data.cygeo
    scale = data.geo_scale
    return plotgeorbf_(shp, coefs, xcoord, ycoord, cx, cy, scale, title, clim)
end

function plotgeorbf_(shp::DataFrame, coefs::Vector{Float64},
                     xcoord::Vector{Float64},
                     ycoord::Vector{Float64},
                     cx::Vector{Float64},
                     cy::Vector{Float64},
                     scale::Float64,
                     title::String,
                     clim)
    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    geo = [interpolant(rbf, xcoord[i], ycoord[i],
                       coefmat(coefs ./ 10, length(cx), length(cy)),
                       cx, cy, scale) for i in eachindex(xcoord)];
    geo = geo .- mean(geo)
    dfgeo = DataFrame(; xcoord, ycoord, geo)
    dfgeo.xcoord = round.(dfgeo.xcoord, digits = 4)
    dfgeo.ycoord = round.(dfgeo.ycoord, digits = 4)
    shp2 = leftjoin(shp, dfgeo, on = [:xcoord, :ycoord])
    shp2 = GeoTable(shp2)
    p = viz(shp2.geometry, color = shp2.geo)
    return dfgeo,  p
end

function plotdensrbf(n::NamedTuple, clim)
    coefs = extract_coefs(n.chn[end, :, 1], "ζ")
    data = n.mdl.mdl.args
    y = n.chn[:year][1]
    return plotdensrbf(coefs, data, "Year $y", clim)
end

function plotdensrbf(coefs::Vector{Float64}, data::NamedTuple,
                     title = nothing, clim = nothing)
    cx = data.cx
    cy = data.cy
    scale = data.rbf_scale
    R = data.R
    return plotdensrbf_(coefs, cx, cy, scale, R, title, clim)
end

function plotdensrbf_(coefs::Vector{Float64},
                      cx::Vector{Float64},
                      cy::Vector{Float64},
                      scale::Float64,
                      R::Vector{Float64}, title = nothing, clim = nothing)
    Rmin, Rmax = extrema(R)
    vals = range(Rmin, Rmax, 1000)
    ## reversing y-values so we don't have to use heatmap(.., yflip =
    ## true), because this also flips the axis labels
    mat = [interpolant(rbf, Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale)
           for Rfrom in vals, Rto in vals]';
    mat = mat .- mean(mat)
    p = heatmap(vals, vals, mat, title = title, clim = clim, aspect_ratio = :equal)
    return mat, p
end
