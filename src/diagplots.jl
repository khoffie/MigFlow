function plotfit(flows, preds)
    idx = subset(flows, 10^3)
    x = log.(preds)[idx]
    y = log.(flows)[idx]
    df = sort(DataFrame(; x, y), :x)
    scatter(df.x, df.y, xlab = L"\log(\hat{y})", ylab = L"\log(y)",
            label = "", alpha = 0.5)
    diagonal!(df.x, df.y)
    smoother!(x, y)
end

function _plotdist(f, flows, preds, dist, lbl)
    idx = subset(flows, 10^3)
    y = log.(flows ./ preds)[idx]
    x = dist[idx]
    df = sort(DataFrame(; x, y), :x)
    f(df.x, df.y,
      xlab = "Distance in km",
      ylab = L"\log(y / \hat{y})", label = lbl, alpha = .5)
    hline!([0], color = :darkred, linewidth = 2, label = "")
    smoother!(x, y)
end

function plotdist(df, preds, lbl = "")
    _plotdist(scatter, df.flows, df[!, preds], df.dist, lbl)
end

function plotdist!(df, preds, lbl = "")
    _plotdist(scatter!, df.flows, df[!, preds], df.dist, lbl)
end

function plotpop(flows, preds, frompop, topop)
    idx = subset(flows, 10^3)
    x = (log.(frompop) .+ log.(topop))[idx]
    y = (log.(flows ./ preds))[idx]
    df = sort(DataFrame(; x, y), :x)
    scatter(df.x, df.y,
         xlab = L"\log(Pop_d \cdot Pop_o)", label = "",
         ylab = L"\log(y / \hat{y})", alpha = .5)
    hline!([0], color = :darkred, linewidth = 2, label = "")
    smoother!(x, y)
end

function plotdens(flows, preds, fromdens, todens)
    flows = DataFrame(flows = flows, preds = preds,
                      fromdens = fromdens, todens = todens)
    ## prevent overplotting
    nrow(flows) > 1e3 ? p = 1000 / nrow(flows) : p = 1
    flows = sample_rows(flows, p)

    plot(flows.fromdens .+ flows.todens,
         log.(flows.flows ./ flows.preds), seriestype = :scatter,
         xlab = L"fromlrd * tolrd", label = "",
         ylab = L"\log(y / \hat{y})")
end

function plotnet(net)
    df = sort(DataFrame(nmrp = net.nmrp, nmr = net.nmr), :nmrp)
    scatter(df.nmrp, df.nmr,
            xlab = "netpred / (influxpred + outfluxpred)",
            ylab = "net / (influx + outflux)",
            xlim = (-1, 1),
            ylim = (-1, 1),
            alpha = .5, label = "")
    diagonal!(df.nmrp, df.nmr)
    smoother!(df.nmrp, df.nmr)
end

function resid(flows, f = mean)
    rs = ((flows.flows - flows.preds).^2) ./ flows.preds
    r = round(f(rs), digits = 0)
    return r
end

function distance_kde(x::AbstractVector, w::AbstractVector,
                      kernel, lbl = nothing, add = false)
    xs = StatsBase.wsample(x, Weights(w), 10^4)
    observed = kde(xs, kernel)
    f = add == true ? plot! : plot
    p = f(observed, xlab = "Distance", ylab = "Density",
##             title = "Migration Distances",
             lw = 2, label = lbl,
          xrange = (minimum(x) - 10, maximum(x) + 10))
    return p
end

subset(x, n) = StatsBase.sample(1:length(x), n)

function pma!(x, y, w)
    plot!(x, moving_average(y, w), linewidth = 5,
      colour = :blue, label = "SMA, m = $w")
end

function diagonal!(x, y)
    xmin, xmax = mm(x)
    ymin, ymax = mm(y)
    low = min(xmin, ymin)
    high = max(xmax, ymax)

    plot!([low, high], [low, high],
          color = :darkred,
          linewidth = 2,
          label = "")
end

function smoother!(x, y, span = .5)
    us = range(extrema(x)...; step = .1)
    vs = Loess.predict(loess(x, y; span = span), us)
    return plot!(us, vs, legend = false, color = "red", linewidth = 4)
end
