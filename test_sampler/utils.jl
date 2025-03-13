logistic(x) = 1 / (1 + exp(-x))

function wssr(flows, preds, f = mean)
    rs = ((flows.flows - flows[:, preds]).^2) ./ flows[:, preds]
    r = round(f(rs), digits = 0)
    return r
end
