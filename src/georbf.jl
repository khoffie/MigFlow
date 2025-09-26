function plotgeo(r::EstimationResult, shp::GeoTable, st::GeoTable)
    geodf = getgeo(r)
    a, y = getageyear(r)
    fig = Figure(size = (400, 400), fontsize = 10)
    yearmap(geodf, shp, st, a, y, fig)
    prettytitle!(fig, "Locational Asymmetries, $a")
    Colorbar(fig[end + 1, :], colorrange = extrema(geodf.geo), vertical = false, height = 3,
             width = Relative(.5))
    return geodf, fig
end

function overlay_states(ax, stshp)
    viz!(ax, stshp.geometry, showsegments = true, alpha = 0,
         segmentcolor = :white, segmentsize = 0.5)
    return ax
end

function getgeo(r::EstimationResult)
    coefs = extract_coefs(r.chn[end, :, 1], "η")
    data = r.mdl.mdl.args
    a, y = getageyear(r)

    xcoord = data.xcoord
    ycoord = data.ycoord
    cx = data.cxgeo
    cy = data.cygeo
    scale = data.geo_scale

    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    geo = [interp(xcoord[i], ycoord[i],
                  coefmat(coefs ./ 10, length(cx), length(cy)),
                  cx, cy, scale) for i in eachindex(xcoord)];
    geo = geo .- mean(geo)
    dfgeo = DataFrame(; xcoord, ycoord, geo)
    dfgeo.agegroup .= a
    dfgeo.year .= y
    return dfgeo[!, [:agegroup, :year, :geo, :xcoord, :ycoord]]
end

function joingeometry(geo, shp)
    shp2 = shp[!, [:x, :y, :geometry]]
    shp2.xcoord = round.(scale_to_unit(shp2.x), digits = 4)
    shp2.ycoord = round.(scale_to_unit(shp2.y), digits = 4)
    geo2 = copy(geo)
    geo2.xcoord = round.(geo2.xcoord, digits = 4)
    geo2.ycoord = round.(geo2.ycoord, digits = 4)
    leftjoin!(geo2, shp2[!, [:xcoord, :ycoord, :geometry]],
              on = [:xcoord, :ycoord])
    return GeoTable(geo2)
end
