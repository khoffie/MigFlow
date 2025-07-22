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

function bound(a, ndc, ngcx, ngcy)
    lbdensity = fill(-100.0, ndc)
    lbgeo = fill(-100.0, ngcx * ngcy)
    ubdensity = fill(100.0, ndc)
    ubgeo = fill(100.0, ngcx * ngcy)

    pastelb(c) = vcat(c, lbdensity..., lbgeo...)
    pasteub(c) = vcat(c, ubdensity..., ubgeo...)
    ## ub alpha only makes sense for distscale = 100 and pop /
    ## median(pop). Otherwise base prob to migrate might be very different
    a == "below18" && return pastelb([-10.0, 10.0, 1.0]), pasteub([-5.0, 40.0, 50.0])
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
