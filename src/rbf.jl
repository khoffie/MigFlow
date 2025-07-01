# Define the RBF function, valid only for r < 1

rbf(x) = abs(x) < one(x) ? exp(one(x) - one(x) / (one(x) - x^2)) : zero(x)
rbfscale(cx, cy, k) = k * LinearAlgebra.norm([cx[1],cy[1]] - [cx[2],cy[2]])
function interpolant(f, x, y,
                     w,
                     cx,
                     cy,
                     scale,
                     N = length(cx))
    @assert length(w) == length(cx) * length(cy)
    res = zero(x)
    for i in eachindex(cx)
        for j in eachindex(cy)
            r = sqrt((x - cx[i])^2 + (y - cy[j])^2) / scale
            res += w[(j - 1) * N + i] * f(r)
        end
    end
    return res
end

function scale_to_unit(x)
    xmin, xmax = extrema(x)
    return [2 * (xi - xmin) / (xmax - xmin) - 1  for xi in x]
end

# districts = year(CSV.read("../data/districts.csv", DataFrame), 2017)

# R = scale_to_unit(log.(districts.density ./ median(districts.density)))
# Rmin, Rmax = extrema(R)
# vals = range(Rmin, Rmax, 1000)
# s = 5
# cx = [range(Rmin, Rmax, s);]
# cy = [range(Rmin, Rmax, s);]
# scale = rbfscale(cx, cy, 1.0)
# w = zeros(Float64, s, s)
# w = zeros(Float64, s^2)
# w[3] = 1

# mat = [interpolant(rbf, xi, yi, w, cx, cy, scale) for xi in vals, yi in vals];
# heatmap(mat)
