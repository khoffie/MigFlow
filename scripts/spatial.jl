using GeoStats, GeoIO, CairoMakie, CSV, DataFrames, Serialization, Turing
using Revise

include("../src/rbf.jl")
include("../src/utils.jl")
include("../src/analysis.jl")
includet("../src/diagrbf.jl")
include("../src/diagchain.jl")
include("../src/model_helpers.jl")
include("../src/estimation.jl")

shp = DataFrame(GeoIO.load("../data/clean/shapes/districts_ext.shp"))
data = gendata([Symbol.("agebelow18")]);
dfgeo = reduce(vcat, loopstruct(data, getgeo))

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

years = [2000, 2001, 2002, 2015, 2016, 2017]
clim = extrema(filter(row -> row.year ∈ years, dfgeo).geo)
years = reshape(years, (3, 2))'

fig = Figure()
for i in 1:size(years, 1)  # row index (subplot row)
    for j in 1:size(years, 2)  # column index (subplot column)
        yr = years[i, j]
        g, fig = plotgeo(year(age(dfgeo, "below18"), yr), shp, fig, i, j, clim)
    end
end

Colorbar(fig[1:2, 4], limits = clim, vertical = true)

display(fig)
