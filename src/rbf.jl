# Define the RBF function, valid only for r < 1

rbf(x) = abs(x) < one(x) ? exp(one(x) - one(x) / (one(x) - x^2)) : zero(x)
rbfscale(cx, cy, k) = k * LinearAlgebra.norm([cx[1],cy[1]] - [cx[2],cy[2]])
function interpolant(f, x, y,
                     w,
                     cx,
                     cy,
                     scale)
    @assert length(w) == length(cx) * length(cy)
    res = zero(x)
    for i in eachindex(cx)
        for j in eachindex(cy)
            r = sqrt((x - cx[i])^2 + (y - cy[j])^2) / scale
            res += w[(j - 1) * length(cx) + i] * f(r)
        end
    end
    return res
end

function scale_to_unit(x)
    xmin, xmax = extrema(x)
    return [2 * (xi - xmin) / (xmax - xmin) - 1  for xi in x]
end

# x = rand(10)
# y = rand(10)
# w = rand(9)
# cx = [.1, .2, .3]
# cy = [.1, .2, .3]
# interpolant(rbf, x[1], y[1], w, cx, cy)
