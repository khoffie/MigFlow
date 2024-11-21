function densitychains(chain::Chains, flows::DataFrame, districts::DataFrame, n::Int,
                       dens_th, title = nothing)
    fromdens, todens, densmin, densmax = densodensd(flows, districts, n)
    plts = [densheatmap(chain[:, :, i], fromdens, todens, densmin, densmax, dens_th)[1]
            for i in 1 : size(chain)[3]]
    p = Plots.plot(plts..., plot_title = title == nothing ? "" : title)
    return p
end

function densheatmap(chain, fromdens, todens, densmin, densmax, dens_th)
    pars = OrderedDict(zip(string.(chain.value.axes[2]), chain.value.data[end, :, 1]))
    kds = [k for k in keys(pars) if contains(k, "kd")]
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
                getindex.(Ref(pars), kds) ./ 10)
    ikd = InterpKDE(kde((fromdens, todens)))
    interval = densmin : .1 : densmax
    p = heatmap(interval, interval, [kdfun(x, y) * (pdf(ikd, x, y) > dens_th)
                                     for y in interval, x in interval])
    return p, values
end

function densodensd(flows, districts, n)
    densities = districts[:, [:distcode, :density]]
    densities.density = log.(districts.density / median(districts.density))
    minlrd = minimum(densities.density)
    maxlrd = maximum(densities.density)
    df = innerjoin(flows, densities, on = :fromdist => :distcode)
    rename!(df, :density => "fromdens")
    leftjoin!(df, densities, on = :todist => :distcode)
    rename!(df, :density => "todens")

    df1 = df[shuffle(1 : nrow(df))[1 : n], :]
    df1.todens = collect(skipmissing(df1.todens)) ## Vector{Union{Missing, Float64}}
    return df1.fromdens, df1.todens, minlrd, maxlrd
end

