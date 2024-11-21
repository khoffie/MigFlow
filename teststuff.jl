using Revise
includet("main.jl")
using KernelDensity

path = "./manuscript_input/30kMH"
f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
districts, flows = loadallGermData(sample = false, positive_only = true)
flows = filter(:agegroup => ==("30-50"), flows)

p = densitychains(chain, flows, districts, 1000)
p

lrd = sort(log.(districts.density / median(districts.density)))
minlrd = minimum(lrd)
maxlrd = maximum(lrd)


pars = OrderedDict(zip(string.(chain.value.axes[2]), chain.value.data[end, :, 1]))
kds = [k for k in keys(pars) if contains(k, "kd")]
kdfun = Fun(ApproxFun.Chebyshev(minlrd .. maxlrd) * ApproxFun.Chebyshev(minlrd .. maxlrd),
            getindex.(Ref(pars), kds) ./ 10)
xes = StatsBase.sample(lrd, 100)
yes = StatsBase.sample(lrd, 100)
ikd = InterpKDE(kde((xes, yes)))
interval = minlrd : .1 : maxlrd
p = Plots.heatmap([kdfun * (pdf(ikd, x, y) > 0.1) for x in interval, y in interval])
p = Plots.heatmap([1 * (pdf(ikd, x, y) > 0.1) for x in interval, y in interval])

param(chain, :lp, 50)
