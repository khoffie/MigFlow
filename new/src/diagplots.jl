function plotfit(flows, preds)
    flows = DataFrame(flows = flows, preds = preds)
    ## prevent overplotting
    nrow(flows) > 1e3 ? p = 1000 / nrow(flows) : p = 1
    flows = sample_rows(flows, p)

    r = resid(flows, mean)

    Plots.plot(log.(flows.preds), log.(flows.flows),
         seriestype = :scatter,
         xlab = L"\log(\hat{y})", ylab = L"\log(y)",
         label = "")
end

function plotdist(flows, preds, dist)
    flows = DataFrame(flows = flows, preds = preds, dist = dist)
    ## prevent overplotting
    nrow(flows) > 1e3 ? p = 1000 / nrow(flows) : p = 1
    flows = sample_rows(flows, p)

    y = log.(flows.flows ./ flows.preds)
    x = flows.dist
    plot(x, y, seriestype = :scatter,
         xlab = L"km",
         ylab = L"\log(y / \hat{y})", label = "")
end

function plotpop(flows, preds, frompop, topop)
    flows = DataFrame(flows = flows, preds = preds,
                      frompop = frompop, topop = topop)
    ## prevent overplotting
    nrow(flows) > 1e3 ? p = 1000 / nrow(flows) : p = 1
    flows = sample_rows(flows, p)

    plot(log.(flows.frompop) .+ log.(flows.topop),
         log.(flows.flows ./ flows.preds), seriestype = :scatter,
         xlab = L"\log(Pop_d \cdot Pop_o)", label = "",
         ylab = L"\log(y / \hat{y})")
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
    scatter(net.nmrp, net.nmr,
            xlab = "netpred / (influxpred + outfluxpred)",
            ylab = "net / (influx + outflux)",
            xlim = (-1, 1),
            ylim = (-1, 1))
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
