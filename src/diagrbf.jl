function plotgeorbf(coefs::Vector{Float64},
                    data::NamedTuple)
    xcoord = data.xcoord
    ycoord = data.ycoord
    cx = data.cxgeo
    cy = data.cygeo
    k = data.kgeo
    return plotgeorbf_(coefs, xcoord, ycoord, cx, cy, k)
end

function plotgeorbf_(coefs::Vector{Float64},
                    xcoord::Vector{Float64},
                    ycoord::Vector{Float64},
                    cx::Vector{Float64},
                    cy::Vector{Float64},
                    k::Float64)
    xmin, xmax = extrema(xcoord)
    ymin, ymax = extrema(ycoord)
    geoscale   = rbfscale(cx, cy, k)
    geo = [interpolant(rbf, xcoord[i], ycoord[i], coefmat(coefs ./ 10),
                       cx, cy, geoscale) for i in eachindex(xcoord)];
    geo = geo .- mean(geo)
    dfgeo = DataFrame(; xcoord, ycoord, geo)
    ratio = (ymax - ymin) / (xmax - xmin)
    ratio = 1.33
    width = 600
    p = scatter(xcoord, ycoord,
                marker_z = dfgeo.geo,
                markersize = 10, size = (600, width * ratio),
                label = "", yflip = false)
    res = collect(IterTools.product(cx, cy))
    xvals = getindex.(res, 1)
    yvals = getindex.(res, 2)
    scatter!(vec(xvals), vec(yvals), color = :red)
    return dfgeo, p
end

function plotdensrbf(coefs::Vector{Float64}, data::NamedTuple)
    cx = data.cx
    cy = data.cy
    k = data.kdens
    R = data.R
    return plotdensrbf_(coefs, cx, cy, k, R)
end

function plotdensrbf_(coefs::Vector{Float64},
                      cx::Vector{Float64},
                      cy::Vector{Float64},
                      k::Float64,
                      R::Vector{Float64})
    Rmin, Rmax = extrema(R)
    vals = range(Rmin, Rmax, 1000)
    s = Int(sqrt(length(coefs)))
    scale = rbfscale(cx, cy, k)
    mat = [interpolant(rbf, Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale)
           for Rto in vals, Rfrom in vals];
    mat = mat .- mean(mat)
    p = heatmap(mat)
    return mat, p
end
