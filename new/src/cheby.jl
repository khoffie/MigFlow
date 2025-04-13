function chebynodes(n, a = -1, b = 1)
    k = 1:(n + 1)                          # Indices
    d = (2k .- 1) ./ (2(n + 1))             # Fractions of π
    x = -cos.(d .* π)                       # Chebyshev nodes in [-1,1]
    x = (a + b) / 2 .+ ((b - a) / 2) .* x    # Scale to [a,b]
    return DataFrame(; k, d, x)
end

function plotgermancheby(pshp, n, districts)
    xmin, xmax = mm(districts.xcoord)
    ymin, ymax = mm(districts.ycoord)
    xnodes = chebynodes(sqrt(n), xmin, xmax)
    ynodes = chebynodes(sqrt(n), ymin, ymax)
    m = collect(Iterators.product(xnodes.x, ynodes.x))
    xvals = [p[1] for p in m]
    yvals = [p[2] for p in m]
    w = 600
    r = (ymax - ymin) / (xmax - xmin)
    line!(x, y) = plot!(x, y, linewidth = 2, color = "blue", label = "")
    p = scatter!(pshp, xvals, yvals, label = "", size = (w, w *r),
            color = "black",
            xlim = (xmin - 10, xmax + 10),
            ylim = (ymin - 10, ymax + 10),
            title = "$n Cheby nodes")
    line!([xmin, xmax], [ymin, ymin])
    line!([xmin, xmax], [ymax, ymax])
    line!([xmin, xmin], [ymax, ymin])
    line!([xmax, xmax], [ymin, ymax])
    return p
end
