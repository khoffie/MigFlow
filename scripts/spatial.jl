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

years = [2006, 2007, 2008, 2015, 2016, 2017]
plotgeoyears(dfgeo, shp, "below18", years)
