# Define the RBF function, valid only for r < 1

rbf(x) = abs(x) < 1 ? exp(1 - 1 / (1 - x^2)) : 0.0

function interpolant(f, x::Float64, y::Float64,
                     w,
                     cx::Vector{Float64},
                     cy::Vector{Float64},
                     k::Int = 2)
    @assert length(w) == length(cx) * length(cy)
    res = 0.0
    for i in eachindex(cx)
        for j in eachindex(cy)
            r = sqrt((x - cx[i])^2 + (y - cy[j])^2) / k
            res += w[i + j] * r
        end
    end
    return res
end

function scale_to_unit(x)
    xmin, xmax = extrema(x)
    return [2 * (xi - xmin) / (xmax - xmin) - 1  for xi in x]
end

x = rand(10)
y = rand(10)
w = rand(9)
cx = [.1, .2, .3]
cy = [.1, .2, .3]
interpolant(rbf, x[1], y[1], w, cx, cy)
