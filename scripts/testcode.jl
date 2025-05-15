using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings
using BenchmarkTools, ADTypes, ReverseDiff, PrettyTables, GLM

include("../src/utils.jl")
include("../src/othermodels.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/fullmodel.jl")
include("../src/fullmodel2.jl")
include("../src/chebies.jl")
include("../src/choiceset.jl")
include("../src/gen_mdat.jl")

data = load_data("30-50", 2017, 0.1, "../data/");

df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
di = addlrd!(CSV.read("../data/districts.csv", DataFrame))
df = joinlrd(df, di)
data = (df = age(year(df, 2017), "30-50"), districts = year(di, 2017))


mdat = gen_mdat(data; type = "joint", distscale = 100.0, ndc = 28, ngc = 1);
out1 = @time estimate(norm, mdat);
out2 = @time estimate(choice, mdat);
out3 = @time estimate(distonly, mdat);
out4 = @time estimate(choice, mdat; normalize = false);
out1.out
out2.out
out3.out
out4.out

mdat = gen_mdat(data; type = "conditional", distscale = 100.0, ndc = 28, ngc = 1);
out1 = @time estimate(norm, mdat);
out2 = @time estimate(choice, mdat);
out3 = @time estimate(distonly, mdat);
out4 = @time estimate(choice, mdat; normalize = false);
out1.out
out2.out
out3.out
out4.out

out.plt[5]

df = out.dens
df = df[df.fromdens .< 2.7 .&& df.todens .< 2.7, :]
heatmap(reshape(df.funval, (90, 90)))

?heatmap

plot(Gamma(1, 1))
names(out1.net)



lcds = unique(sort(DataFrame(; mdat.fromdist, mdat.fromdist_orig), :fromdist))
leftjoin!(out1.net, lcds, on = :fromdist)
rename!(out1.net, :fromdist_orig => :distcode)

sort(out1.net, :nmrp)[!, [:distcode, :influx, :influxp, :outflux, :outfluxp, :nmr, :nmrp]]
scatter(out1.net.nmrp, out1.net.nmr)
