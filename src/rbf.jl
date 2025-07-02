# Define the RBF function, valid only for r < 1

rbf(x) = abs(x) < one(x) ? exp(one(x) - one(x) / (one(x) - x^2)) : zero(x)
rbfscale(cx, cy, k) = k * LinearAlgebra.norm([cx[1],cy[1]] - [cx[2],cy[2]])

function interpolant(f, x, y, w, cx, cy, scale)
    res = zero(x)
    @inbounds for i in eachindex(cx), j in eachindex(cy)
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

xmin = 0
xmax = 1
x = range(xmin, xmax, 100)
y = x
N = 81
s = Int(sqrt(N))
w = coefmat(randn(N))
w = zeros((s, s))

w[5, 5] = 1
cx = [range(0, 1, s);]
cy = [range(0, 1, s);]
rbf_scale = rbfscale(cx, cy, 1.0)
interpolant(rbf, x[1], y[1], w, cx, cy, rbf_scale)
@benchmark interpolant(rbf, x[1], y[1], w, cx, cy, rbf_scale)

mat = [interpolant(rbf, xi, yi, w, cx, cy, rbf_scale) for xi in x, yi in y];
heatmap(mat .- mean(mat))
