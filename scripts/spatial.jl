using GeoStats, GeoIO, CairoMakie, CSV, DataFrames, Serialization, Turing
using Revise

include("../src/rbf.jl")
include("../src/utils.jl")
include("../src/analysis.jl")
includet("../src/diagrbf.jl")
include("../src/diagchain.jl")
include("../src/model_helpers.jl")
include("../src/estimation.jl")

shp = GeoIO.load("../data/clean/shapes/districts_ext.shp")
st = GeoIO.load("../data/clean/shapes/states.shp")
typeof(st)
fig = Figure()
plotgeo_(year(dfgeo, 2001), DataFrame(shp), fig, nothing, nothing, nothing)

years = [2006, 2007, 2008, 2015, 2016, 2017]
F1 = plotgeoyears(dfgeo, DataFrame(shp), st, "below18", years)
display(F1)
