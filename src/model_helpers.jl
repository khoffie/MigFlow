rbf(x) = abs(x) < one(x) ? exp(one(x) - one(x) / (one(x) - x^2)) : zero(x)
rbfscale(cx, cy, k) = k * LinearAlgebra.norm([cx[1],cy[1]] - [cx[2],cy[2]])

function interpolant(f, x, y, w, cx, cy, scale)
    res = zero(x)
    @inbounds for i in eachindex(cx), j in reverse(eachindex(cy))
        dx = x - cx[i]
        dy = y - cy[length(cy) + 1 - j]
##         dy = y - cy[j]
        r = sqrt(dx * dx + dy * dy) / scale
        res += w[j, i] * f(r)
    end
    return res
end

# ### Check density RBF
# vals = range(-1, 1, 1000)
# cx = collect(range(-1, 1, 4))
# cy = collect(range(-1, 1, 4))

# w = [0, 1, 0, 0,
#      0, 0, 0, 0,
#      0, 0, 2, 0,
#      1, 0, 0, 0]

# extract_coefs(out.chn, "ζ")

# coefmat(w)
# mat, p = plotdensrbf_(Float64.(w), cx, cy, .5, Float64.(vals), "test", nothing);

function scale_to_unit(x)
    xmin, xmax = extrema(x)
    return [2 * (xi - xmin) / (xmax - xmin) - 1  for xi in x]
end

fdist(D, ds) = D / ds
fradius(P, ρ, ds) = sqrt((P / ρ) / 2π) / ds
lc(x) = levelcode.(categorical(x))
StMvN(n, σ) = MvNormal(zeros(n), fill(σ, n))

function coefmat(c)
    s = Int(sqrt(length(c)))
    return reshape(c, (s, s))'
end

function coefmat(c, nx::Int, ny::Int)
    return reshape(c, (nx, ny))'
end

function genfrompop(df, type)
    type == "joint" && return df.frompop
    df2 = combine(groupby(data.df, :fromdist), :flows => sum)
    return leftjoin(df, df2, on = :fromdist).flows_sum
end

function geo_ratio(xcoord, ycoord, xlim, ylim)
    ## For equidistant spacing of Geo RBF centers, we need more
    ## coordinates in y direction, because Germany is about 1.33 times
    ## longer than wide.
    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    ratio = (ymax - ymin) / (xmax - xmin)
    ## Also we do not need all the grid, so we want to adjust for that
    used_x = (xlim[2] - xlim[1]) / 2
    used_y = (ylim[2] - ylim[1]) / 2
    return ratio * (used_y / used_x)
end

function bound(a, ndc, ngcx, ngcy)
    lbdensity = fill(-100.0, ndc)
    lbgeo = fill(-100.0, ngcx * ngcy)
    ubdensity = fill(100.0, ndc)
    ubgeo = fill(100.0, ngcx * ngcy)

    pastelb(c) = vcat(c, lbdensity..., lbgeo...)
    pasteub(c) = vcat(c, ubdensity..., ubgeo...)
    ## ub alpha only makes sense for distscale = 100 and pop /
    ## median(pop). Otherwise base prob to migrate might be very different
    a == "below18" && return pastelb([-12.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 50.0])
    a == "18-25" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 30.0])
    a == "25-30" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 30.0])
    a == "30-50" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
    a == "50-65" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
    a == "above65" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 40.0])
end

rbfinits(N, σ, t = 90.0) = clamp.(rand(MvNormal(zeros(N), σ^2 *I(N))), -t, t)
function initialize(a, ndc, ngcx, ngcy)
    density = rbfinits(ndc, 40.0)
    geo = rbfinits(ngcx * ngcy, 10.0)
    a == "below18" && return [-8.0, 25.0, 30.0, density..., geo...]
    a == "18-25" && return [-7.0, 18.0, 20.0, density..., geo...]
    a == "25-30" && return [-7.0, 18.0, 20.0, density..., geo...]
    a == "30-50" && return [-6.5, 20.0, 18.0, density..., geo...]
    a == "50-65" && return [-7.5, 23.0, 30.0, density..., geo...]
    a == "above65" && return [-7.5, 23.0, 30.0, density..., geo...]
end
