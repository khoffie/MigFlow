using DataFrames, CSV, CairoMakie, Statistics, StatsBase, Revise
using LaTeXStrings, AlgebraOfGraphics, GeoIO, GeoStats, GLM
using KernelDensity
outp = "/home/konstantin/paper/sections/texdocs/images/"

includet("../src/diagplots.jl")
includet("../src/analyzegeo.jl")
includet("../src/model_helpers.jl")
includet("../src/samplecircle.jl")

df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
df = rmreg(df, :fromdist)
df = rmreg(df, :todist)

di = CSV.read("../data/districts.csv", DataFrame)
di.area = di.pop ./ di.density;

pop = combine(groupby(unique(df, [:fromdist, :year, :agegroup]),
                      [:year, :agegroup]), :frompop => sum => :agepop)
## combine(groupby(pop, :year), :agepop => sum)
shp = GeoIO.load("../data/clean/shapes/districts_ext.shp");
st = GeoIO.load("../data/clean/shapes/states.shp")

years = unique(df.year)
ages = unique(df.agegroup)


######################################################################
################# total flows per region  ############################
net = calc_net(df, :flows);
net = rmreg(net, :fromdist)
net.geo = log.(net.total ./ length(years));
net.yearly_total = net.total ./ length(years)
net.yearly_influx = net.influx ./ length(years)
net.yearly_outflux = net.outflux ./ length(years)

diy = combine(groupby(di, :distcode), [
    ## should pop be pop of Germans only?
    :pop => mean => :pop,
    :area => first => :area,
    :name => first => :name
])

leftjoin!(net, diy, on = :fromdist => :distcode)
net.flow_share = net.yearly_total ./ net.pop .* 100
net.influx_share = net.yearly_influx ./ net.pop .* 100
net.outflux_share = net.yearly_outflux ./ net.pop .* 100

## sort(net, :flow_share)[!, [:name, :flow_share]][350:400, :]
## Makie.density(net.flow_share)
# Makie.scatter(net.area, net.influx_share)
# Makie.scatter(net.area, net.outflux_share)

fig = Figure(size = (400, 400), fontsize = 10);
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title = "Influx + Outflux",
          subtitle = "All age group, yearly average 2000 - 2017")
viz!(ax, rmreg(shp, :AGS).geometry, color = net.geo, colorrange = extrema(net.geo));
overlay_states(ax, st)

ax2 = Axis(fig[1, 2], aspect=DataAspect(),
          title = "(Influx + Outflux) / Population",
          subtitle = "All age group, yearly average 2000 - 2017")
viz!(ax2, rmreg(shp, :AGS).geometry, color = net.flow_share, colorrange = extrema(net.flow_share));
overlay_states(ax2, st)

Colorbar(fig[2, 1], limits = extrema(net.geo), vertical = false,
         height = 3, width = Relative(.4), label = "log(influx + outflux)")
Colorbar(fig[2, 2], limits = extrema(net.flow_share), vertical = false,
         height = 3, width = Relative(.4), label = "(influx + outflux) / pop")

hideall!(ax)
hideall!(ax2)
fig
resize_to_layout!(fig)
save(joinpath(outp, "totalmap.pdf"), fig)

# fig = Figure(size = (400, 400), fontsize = 10);
# ax = Axis(fig[1, 1], aspect=DataAspect(),
#           title = "Influx / Pop",
#           subtitle = "All age group, yearly average 2000 - 2017")
# viz!(ax, rmreg(shp, :AGS).geometry, color = net.influx_share, colorrange = extrema(net.influx_share));
# overlay_states(ax, st)

# ax2 = Axis(fig[1, 2], aspect=DataAspect(),
#           title = "Outflux / Population",
#           subtitle = "All age group, yearly average 2000 - 2017")
# viz!(ax2, rmreg(shp, :AGS).geometry, color = net.outflux_share, colorrange = extrema(net.outflux_share));
# overlay_states(ax2, st)

# Colorbar(fig[2, 1], limits = extrema(net.influx_share), vertical = false,
#          height = 3, width = Relative(.4), label = "influx / pop")
# Colorbar(fig[2, 2], limits = extrema(net.outflux_share), vertical = false,
#          height = 3, width = Relative(.4), label = "outflux / pop")

# hideall!(ax)
# hideall!(ax2)
# fig
