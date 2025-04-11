function t(x, c, m = 100.0, e = 0.1)
    B = (c + 1) / ((m + e)^(c+1) - e^(c+1))
    return B * (x + e)^c
end

## quadgk(x -> t(x, 3, m, 0.1), 0, m, rtol=1e-8)

area(x) = pi * x^2

function cdfdf(x, g, c)
    gc(x) = g(x, c)
    y = gc.(x)
    int = [quadgk(x -> gc(x), 0, x[i], rtol=1e-8) for i in eachindex(x)]
    cdf = first.(int)
    a = area.(x)
    return DataFrame(; x, a, y, cdf)
end

function plotcdf(x, g, c, plt)
    df = cdfdf(x, g, c)
    lbl = "c = $c"

    p1, p2, p3, p4 = plt
    p1 = plot!(p1, df.x, df.y, xlab = "distance, km", ylab = "pdf", label = "")
    p2 = plot!(p2, df.x, df.cdf, xlab = "distance, km", ylab = "cdf", label = "")
    p3 = plot!(p3, df.a, df.y, xlab = "circle area", ylab = "pdf", label = "")
    p4 = plot!(p4, df.a, df.cdf, xlab = "circle area", ylab = "cdf", label = lbl)
    return df, (p1, p2, p3, p4)
end

function cdfplots()
    m = 25
    x = range(0, m, 100)
    g(x, c) = t(x, c, m, 1.0)
    p = plot(), plot(), plot(), plot()

    df, p1 = plotcdf(x, g, -2.0, p)
    df, p2 = plotcdf(x, g, 2.0, p)
    df, p3 = plotcdf(x, g, 3.0, p)
    display(plot(p3...))
end

