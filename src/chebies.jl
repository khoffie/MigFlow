function defdensitycheby(coefs, densmin, densmax)
    cheby = ApproxFun.Chebyshev(densmin .. densmax)
    return Fun(cheby * cheby, coefs)
end

function evaldens(out, d)
    idx = occursin.("ζ", names(out)[1])
    if sum(idx) > 0
        ## Float, otherwise ./ 100 fails, even with out[idx].array ./
        ## 100 Appending 0.0, because, the first coef is no parameter,
        ## but fixed. Because it is no parameter it is not returned by the model
        coefs = vcat(0.0, Float64.(collect(out[idx])) ./ 100)
        return evaldensitycheby(coefs, d.Rmin, d.Rmax)
    else
        return nothing, nothing
    end
end

function evaldensitycheby(coefs, densmin, densmax,
                          vals = range(densmin, densmax, 100),
                          show_plt = false)
    densitycheby = defdensitycheby(coefs, densmin, densmax)
    df = ((fromdens = fd,
          todens = td,
          funval = densitycheby(fd, td)) for fd in vals, td in vals)
    df = DataFrame(df)
    ## otherwise upper right corner too bright
    df = df[df.fromdens .< 2.8 .&& df.todens .< 2.8, :]
    s = Int(sqrt(nrow(df)))
    p = heatmap(vals[1:s], vals[1:s], reshape(df.funval, (s, s)))
    if show_plt; display(p); end
    return df, p
end

function defgeocheby(coefs, xmin, xmax, ymin, ymax)
    chebyx = ApproxFun.Chebyshev(xmin .. xmax)
    chebyy = ApproxFun.Chebyshev(ymin .. ymax)
    return Fun(chebyx * chebyy, coefs)
end

function evalgeo(out, d)
    idx = occursin.("η", names(out)[1])
    if sum(idx) > 0
        ## 100 Appending 0.0, because, the first coef is no parameter,
        ## but fixed. Because it is no parameter it is not returned by the model
        coefs = vcat(0.0, Float64.(collect(out[idx])) ./ 100)
        return evalgeocheby(coefs, d.distcode, d.xcoord, d.ycoord)
    else
        return nothing, nothing
    end
end

function evalgeocheby(coefs, distcode, xcoord, ycoord,
                      show_plt = false)
    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    geocheby = defgeocheby(coefs, xmin, xmax, ymin, ymax)
    geo = geocheby.(xcoord, ycoord)
    geo2 = exp.(geo)
    ##geo = ForwardDiff.value.(geo)
    df = DataFrame(; distcode, xcoord, ycoord, geo, geo2)

    ratio = (ymax - ymin) / (xmax - xmin)
    width = 600
    p = scatter(df.xcoord, df.ycoord,
                marker_z = df.geo,
                markersize = 10, size = (600, width * ratio),
                label = "")
    if show_plt; display(p); end
    return df, p
end

coeforder(m) = [(j, d - j) for d in 0:m for j in 0:d]
coefindicesx(order, m = 5) = findall(t -> t[1] == order, coeforder(m))
coefindicesy(order, m = 5) = findall(t -> t[2] == order, coeforder(m))
