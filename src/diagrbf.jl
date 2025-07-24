function plotgeoyears(dfgeo, shp, agegroup, years)
    df = age(dfgeo, agegroup)
    clim = extrema(filter(row -> row.year ∈ years, df).geo)
    years = reshape(years, (3, 2))'
    fig = Figure()
    for i in 1:size(years, 1)
        for j in 1:size(years, 2)
            yr = years[i, j]
            g, fig = plotgeo(year(df, yr), shp, fig, i, j, clim)
        end
    end
    Colorbar(fig[1:2, 4], limits = clim, vertical = true)
    return fig
end

function plotgeo(geo, shp, fig, x = 1, y = 1, clim = nothing)
    geo2 = joingeometry(geo, shp)
    age = unique(geo.agegroup)[1]
    year = unique(geo.year)[1]
    if isnothing(clim); clim = extrema(geo.geo); end
    ax = Axis(fig[x, y], aspect=DataAspect(), title = "$year")
    hidedecorations!(ax)
    hidespines!(ax)
    viz!(ax, geo2.geometry, color=geo2.geo, colorrange = clim)
    return geo, fig
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
    p = Plots.heatmap(vals, vals, mat, title = title, clim = clim, aspect_ratio = :equal)
    return mat, p
end
