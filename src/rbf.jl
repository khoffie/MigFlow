# Define the RBF function, valid only for r < 1

rbf(x) = abs(x) < 1 ? exp(1 - 1 / (1 - x^2)) : 0.0

function interpolant(f, x, y,
                     w::Vector{Float64},
                     centers::Vector{Tuple{Float64, Float64}},
                     k = 2)
    res = 0.0
    for (i, (cx, cy)) in enumerate(centers)
        r = sqrt((x - cx)^2 + (y - cy)^2) / k
        res += w[i] * f(r)
    end
    return res
end

function scale_to_unit(x)
    xmin, xmax = extrema(x)
    return [2 * (xi - xmin) / (xmax - xmin) - 1  for xi in x]
end
