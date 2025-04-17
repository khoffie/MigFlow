function calc_net(flows, col)
    netf = combine(groupby(flows, :fromdist), col => sum)
    rename!(netf, string(col) * "_sum" => :influx)
    nett = combine(groupby(flows, :todist), col => sum)
    rename!(nett, string(col) * "_sum" => :outflux)
    net = innerjoin(netf, nett, on = [:fromdist => :todist])
    net.net = net.influx .- net.outflux
    net.total = net.influx .+ net.outflux
    net.nmr = net.net ./ net.total
    return net
end

function calc_net_df(flows)
    net = calc_net(flows, :flows)
    netp = calc_net(flows, :preds)

    new = names(netp) .* "p"
    names(netp)
    rename!(netp, names(netp) .=> new)
    net = innerjoin(net, netp, on = :fromdist => :fromdistp)
    return net
end

function moving_average(A::AbstractArray, m::Int)
    ## https://discourse.julialang.org/t/smoothing-noisy-data-using-moving-mean/65329/6
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = m÷2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    return out
end
