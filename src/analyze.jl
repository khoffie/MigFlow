function analyze(r::EstimationResult)
    df = modeldf(r)
    net = netdf(r)
    quick = quickdf(r)
    asym = asymdf(df)
    fig = Figure(size = (1000, 400), fontsize = 20);
    pointsize = 6
    ax1 = Axis(fig[1, 1],
               xlabel = L"\log \hat{y}",
               ylabel = L"\log y",
               title = L"\text{Mean deviance:}%$(quick.deviance[1])",
               aspect = DataAspect(),
               xgridvisible = false, ygridvisible = false)
    plotfit!(ax1, df.flows, df.preds, pointsize)

    ax2 = Axis(fig[1, 2],
               xlabel = L"(\hat{i} - \hat{o}) / (\hat{i} + \hat{o})",
               ylabel = L"(i - o) / (i + o)",
               title = L"\textrm{skillscore} = %$(round(quick.skillscore[1], digits = 2))",
               aspect = DataAspect(),
               xgridvisible = false, ygridvisible = false)
    Makie.ylims!(ax2, (-1, 1))
    Makie.xlims!(ax2, (-1, 1))
    plotasym!(ax2, net, pointsize)

    ax3 = Axis(fig[1, 3],
               xlabel = L"\text{Distance (km)}",
               ylabel = L"\log(y / \hat{y})",
               xgridvisible = false, ygridvisible = false)
    plotdist!(ax3, df.flows, df.preds, df.dist, pointsize)

    ax4 = Axis(fig[1, 4],
               xlabel = L"\log(A_o  P_d)",
               ylabel = L"\log(y / \hat{y})",
               xgridvisible = false, ygridvisible = false)
    plotpop!(ax4, df.flows, df.preds, df.A, df.P, pointsize)
    return (; df, net, quick, asym, fig)
end

function modeldf(r::EstimationResult)
    a, y = getageyear(r)
    data = r.mdl.mdl.args
    df = DataFrame(
        fromdist = data.from,
        todist = data.to,
        flows = data.Y,
        preds = r.prd,
        dist = 100data.D, ## scaling back to original, better grab ds?
        A = data.A,
        P = exp.(data.P[data.to]) ## bec log(P) is saved
    )
    return df
end

function netdf(r::EstimationResult)
        return add_meta(r, calc_net_df(modeldf(r)))
end

function quickdf(r::EstimationResult)
    a, y = getageyear(r)
    m = getmodel(r)
    df = modeldf(r)
    net = netdf(r)
    dev = round(deviance2(df.flows, df.preds), digits = 2)
    err = 100mae(net.asym, net.asymp)
    trivial = 100mae(net.asym, 0)
    skill = skillscore(net.asym, net.asymp)
    quick = DataFrame(model = m, agegroup = a, year = y, deviance = dev,
                      mae = err, mae0 = trivial, skillscore = skill)
    return quick
end

function calc_net(df, col)
    netf = combine(DataFrames.groupby(df, :fromdist), col => sum)
    rename!(netf, string(col) * "_sum" => :outflux)
    nett = combine(DataFrames.groupby(df, :todist), col => sum)
    rename!(nett, string(col) * "_sum" => :influx)
    net = innerjoin(netf, nett, on = [:fromdist => :todist])
    pop = unique(df, :fromdist)[!, [:fromdist, :A]]
    net = innerjoin(net, pop, on = [:fromdist])
    net.net = net.influx .- net.outflux
    net.total = net.influx .+ net.outflux
    net.asym = net.net ./ net.total
    net.nmra = net.net ./ net.A
    return net
end

function calc_net_df(df)
    net = calc_net(df, :flows)
    netp = calc_net(df, :preds)

    new = names(netp) .* "p"
    names(netp)
    rename!(netp, names(netp) .=> new)
    net = innerjoin(net, netp, on = :fromdist => :fromdistp)
    return net
end

function add_meta(r, df)
    a, y = getageyear(r)
    df.agegroup .= a
    df.year .= y
    df.model .= getmodel(r)
    first = ["model","agegroup", "year"]
    last = setdiff(names(df), first)
    select!(df, vcat(first, last))
end

subset(x, n) = StatsBase.sample(1:length(x), n)

function extract_coefs(chn::Chains)
    return vec(chn.value.data)
end

function extract_coefs(chn::Chains, string::String)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end

function extract_sample(chn, type = "best")
    if type == "best"
        m = argmax(chn[:lp].data)
        chn = chn[m[1], :, m[2]]
    end
    return chn
end

function plotfit!(ax, flows, preds, size)
    idx = subset(flows, 10^3)
    x = log.(preds)[idx]
    y = log.(flows)[idx]
    df = DataFrame(; x, y)
    Makie.scatter!(ax, df.x, df.y, alpha = .5, markersize = size)
    diagonal!(ax, df.x, df.y)
    smoother!(ax, df.x, df.y)
end

function plotasym!(ax, net, size)
    Makie.scatter!(ax, net.asymp, net.asym, alpha = .5, markersize = size)
    diagonal!(ax, net.asymp, net.asym)
    smoother!(ax, net.asymp, net.asym)
end

function plotdist!(ax, flows, preds, dist, size)
    idx = subset(flows, 10^3)
    y = log.(flows ./ preds)[idx]
    x = dist[idx]
    df = sort(DataFrame(; x, y), :x)
    Makie.scatter!(ax, df.x, df.y, alpha = .5, markersize = size)
    Makie.hlines!(ax, [0], color = :darkred, linewidth = 2)
    smoother!(ax, df.x, df.y)
end

function plotpop!(ax, flows, preds, frompop, topop, size)
    idx = subset(flows, 10^3)
    x = (log.(frompop) .+ log.(topop))[idx]
    y = (log.(flows ./ preds))[idx]
    df = sort(DataFrame(; x, y), :x)
    Makie.scatter!(ax, x, y, alpha = .5, markersize = size)
    Makie.hlines!(ax, [0], color = :darkred, linewidth = 2)
    smoother!(ax, x, y)
end

function deviance2(y, p)
    loss = zeros(length(y))
    for i in eachindex(y)
        loss[i] = y[i] * log(y[i] / p[i]) - (y[i] - p[i])
    end
    return 2mean(loss)
end

# mse(y, p) = mean((y .- p) .^ 2)
mae(y, p) = mean(abs.(y .- p))
## skillscore(y, p) = 1 - (mse(y, p) / mse(y, 0))
skillscore(y, p) = 1 - (mae(y, p) / mae(y, 0))

function asymdf(df)
    dfod = select(df, :fromdist, :todist, :flows => :outflux,
                  :preds => :outpreds)
    dfdo = select(df, :fromdist => :todist, :todist => :fromdist,
                  :flows => :influx, :preds => :inpreds)
    df1 = leftjoin(dfod, dfdo, on = [:fromdist, :todist])

    df1.total = df1.influx .+ df1.outflux
    df1.totalp = df1.inpreds .+ df1.outpreds
    df1.asym = (df1.influx .- df1.outflux) ./ df1.total
    df1.asymp = (df1.inpreds .- df1.outpreds) ./ df1.totalp
    return dropmissing!(df1)
end
