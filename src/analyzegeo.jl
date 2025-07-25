function plotgeoyears(dfgeo::DataFrame, shp::GeoTable, st::GeoTable,
                      agegroup::String, years::Vector{Int})
    df = age(dfgeo, agegroup)
    clim = extrema(filter(row -> row.year ∈ years, df).geo)
    years = reshape(years, (3, 2))'
    fig = Figure(size = (600, 400), fontsize = 10)
    for i in 1:size(years, 1)
        for j in 1:size(years, 2)
            yr = years[i, j]
            g, fig = plotgeo(year(df, yr), shp, st, fig, i, j, clim)
        end
    end
    Colorbar(fig[0, :], limits = clim, vertical = false,
             height = 5, width = Relative(.5))
    fig[-1, :] = Label(fig, "$agegroup")
    return fig
end

function plotgeo(geo, shp, stshp, fig, x = 1, y = 1, clim = nothing)
    geo2, ax = plotgeo_(geo, shp, fig, x, y, clim)
    overlay_states(ax, stshp)
    return geo2, fig
end

function plotgeo_(geo, shp, fig, x, y, clim)
    geo2 = joingeometry(geo, DataFrame(shp))
    age = unique(geo.agegroup)[1]
    year = unique(geo.year)[1]
    if isnothing(clim); clim = extrema(geo.geo); end
    ax = Axis(fig[x, y], aspect=DataAspect(), title = "$year")
    tightlimits!(ax)
    hidedecorations!(ax)
    hidespines!(ax)
    viz!(ax, geo2.geometry, color=geo2.geo, colorrange = clim)
    return geo2, ax
end

function overlay_states(ax, stshp)
    viz!(ax, stshp.geometry, showsegments = true, alpha = 0,
         segmentcolor = :white, segmentsize = 0.5)
    return ax
end

function getgeo(n::NamedTuple)
    coefs = extract_coefs(n.chn[end, :, 1], "η")
    data = n.mdl.mdl.args
    y = Int(n.chn[:year][1])
    a = recodeage(Int(n.chn[:age][1]))

    xcoord = data.xcoord
    ycoord = data.ycoord
    cx = data.cxgeo
    cy = data.cygeo
    scale = data.geo_scale

    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    geo = [interpolant(rbf, xcoord[i], ycoord[i],
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
