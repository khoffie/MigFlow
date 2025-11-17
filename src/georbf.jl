function plotgeo(r::EstimationResult, shp::GeoTable, st::GeoTable,
                 crange = nothing,
                 fig = Figure(size = (400, 400), fontsize = 10),
                 x = 1, y = 1, legend = true)
    geodf = getgeo(r)
    a, yr = getageyear(r)
    if isnothing(crange); crange = calc_crange(geodf.geo); end
    plotgeo(geodf, shp, st, crange, a, yr, fig, x, y, legend)
    return geodf, fig
end

function plotgeo(df::DataFrame, shp, st, crange, a, yr, fig,
                 x = 1, y = 1, legend = true)
    df = year(age(df, a), yr)
    df = joingeometry(df, DataFrame(shp))
    ax = Axis(fig[x, y], aspect=DataAspect(), title = "$yr")
    viz!(ax, df.geometry, color = df.geo, colorrange = crange,
         colormap = :roma);
    hideall!(ax)
    overlay_states(ax, st)
    if legend
        prettytitle!(fig, "Locational Asymmetries, $a")
        Colorbar(fig[end + 1, :], colorrange = crange, colormap = :roma,
                 vertical = false, height = 3, width = Relative(.5))
    end
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
