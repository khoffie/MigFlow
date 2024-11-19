includet("src/datafitting.jl")

path = "./manuscript_input/tempered"
f = "germchain_2017_30-50.csv"

chain = deserialize(joinpath(path, f))
geodat, agedat = loadallGermData(sample = false, positive_only = true)
agedat = filter(:agegroup => ==("30-50"), agedat)

logreldens = log.(geodat.density / median(geodat.density))
minlrd = minimum(logreldens)
maxlrd = maximum(logreldens)

function densheatmap(chain, densmin, densmax, logreldens, title, type = nothing )
    pars = OrderedDict(zip(string.(chain.value.axes[2]), chain.value.data[end, :, 1]))
    kds = [k for k in keys(pars) if contains(k, "kd")]
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
                getindex.(Ref(pars), kds) ./ 10)
    p = Plots.heatmap(kdfun, title = type)
    return p
end

p1 = densheatmap(chain, minlrd, maxlrd, "logreldens")
p2 = densheatmap(chain, 0, 5000, "untransformed", "value")
plot(p1, p2)
p1


p = densitychains(chain, flows, districts, 100)
