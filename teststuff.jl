using LogDensityProblems: logdensity
using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots, StatsPlots, LaTeXStrings, Revise

path = "./manuscript_input/tempered"
f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
flows = CSV.read("./data/FlowDataGermans.csv", DataFrame)
districts = CSV.read("./data/districts.csv", DataFrame)

function densodensd(flows, districts, n)
    densities = districts[:, [:distcode, :density]]
    leftjoin!(flows, densities, on = :fromdist => :distcode)
    rename!(df, :density => "fromdens")
    leftjoin!(flows, densities, on = :todist => :distcode)
    rename!(df, :density => "todens")

    df1 = df[shuffle(1 : nrow(df))[1 : n], [:fromdens, :todens]]
    fromdens = sort(df1.fromdens)
    todens = sort(df1.todens)
    return fromdens, todens
end

function densitychains(chain, denso, densd, title = nothing)
    plts = [plotdensity(chain[:, :, i], denso, densd)[1] for i in 1 : size(chain)[3]]
    p = Plots.plot(plts..., plot_title = title == nothing ? "" : title)
    return p
end

function plotdensity(chain, denso, densd)
    densmin = 0
    densmax = 5000
    pars = OrderedDict(zip(string.(chain.value.axes[2]), chain.value.data[end, :, 1]))
    kds = [k for k in keys(pars) if contains(k, "kd")]
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
                getindex.(Ref(pars), kds) ./ 10)
    values = [kdfun(x, y) for x in denso, y in densd]
    p = Plots.heatmap(denso, densd, vals, colorbar = false, ticks = false)
    return p, values
end

p = densitychains(chain, fromdens, todens)


p1, v1 = plotdensity(chain[:, :, 1], fromdens, todens)
p2, v2 = plotdensity(chain[:, :, 2], fromdens, todens)

plot(p1, p2)
p
