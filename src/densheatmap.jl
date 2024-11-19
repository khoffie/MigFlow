function densitychains(chain::Chains, flows::DataFrame, districts::DataFrame, n::Int, title = nothing)
    fromdens, todens = densodensd(flows, districts, n)
    plts = [densheatmap(chain[:, :, i], fromdens, todens)[1] for i in 1 : size(chain)[3]]
    p = Plots.plot(plts..., plot_title = title == nothing ? "" : title)
    return p
end

function densheatmap(chain, denso, densd)
    densmin = 0
    densmax = 5000
    pars = OrderedDict(zip(string.(chain.value.axes[2]), chain.value.data[end, :, 1]))
    kds = [k for k in keys(pars) if contains(k, "kd")]
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
                getindex.(Ref(pars), kds) ./ 10)
    # values = [kdfun(x, y) for x in denso, y in densd]
    # p = Plots.heatmap(denso, densd, values, colorbar = false, ticks = false)
    p = Plots.heatmap(kdfun)
    return p, values
end

function densodensd(flows, districts, n)
    densities = districts[:, [:distcode, :density]]
    df = innerjoin(flows, densities, on = :fromdist => :distcode)
    rename!(df, :density => "fromdens")
    leftjoin!(df, densities, on = :todist => :distcode)
    rename!(df, :density => "todens")

    df1 = df[shuffle(1 : nrow(df))[1 : n], [:fromdens, :todens]]
    fromdens = sort(df1.fromdens)
    todens = sort(df1.todens)
    return fromdens, todens
end
