# Define the RBF function, valid only for r < 1

rbf(x) = abs(x) < one(x) ? exp(one(x) - one(x) / (one(x) - x^2)) : zero(x)
rbf2(x) = exp(- (x / .35) ^ 2)
rbfscale(cx, cy, k) = k * LinearAlgebra.norm([cx[1],cy[1]] - [cx[2],cy[2]])

function interpolant(f, x, y, w, cx, cy, scale)
    res = zero(x)
    @inbounds for i in reverse(eachindex(cx)), j in eachindex(cy)
        dx = x - cx[i]
        dy = y - cy[j]
        r = sqrt(dx * dx + dy * dy) / scale
        res += w[i, j] * f(r)
    end
    return res
end

function scale_to_unit(x)
    xmin, xmax = extrema(x)
    return [2 * (xi - xmin) / (xmax - xmin) - 1  for xi in x]
end
