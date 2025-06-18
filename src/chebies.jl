function defdensitycheby(coefs, densmin, densmax)
    cheby = ApproxFun.Chebyshev(densmin .. densmax)
    return Fun(cheby * cheby, coefs)
end

function defgeocheby(coefs, xmin, xmax, ymin, ymax)
    chebyx = ApproxFun.Chebyshev(xmin .. xmax)
    chebyy = ApproxFun.Chebyshev(ymin .. ymax)
    return Fun(chebyx * chebyy, coefs)
end

function evaldens(out, d)
    if any(occursin.("ζ", names(out)[1]))
        kds = out[["ζ_raw[$i]" for i in 1 : d.ndc]]
        return evaldensitycheby(kds, d.Rmin, d.Rmax)
    else
        return nothing, nothing
    end
end

function evalgeo(out, d)
    if any(occursin.("η", names(out)[1]))
        kgs = out[["η_raw[$i]" for i in 1 : d.ngc]]
        return evalgeocheby(kgs, d.distcode, d.xcoord, d.ycoord)
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
    ##geo = ForwardDiff.value.(geo)
    df = DataFrame(; distcode, xcoord, ycoord, geo)

    ratio = (ymax - ymin) / (xmax - xmin)
    width = 600
    p = scatter(df.xcoord, df.ycoord,
                marker_z = df.geo,
                markersize = 10, size = (600, width * ratio),
                label = "")
    if show_plt; display(p); end
    return df, p
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

coeforder(m) = [(j, d - j) for d in 0:m for j in 0:d]
coefindicesx(order, m = 5) = findall(t -> t[1] == order, coeforder(m))
coefindicesy(order, m = 5) = findall(t -> t[2] == order, coeforder(m))
