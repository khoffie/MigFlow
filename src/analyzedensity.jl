function plotdensrbf(n::NamedTuple, clim)
    coefs = extract_coefs(n.chn[end, :, 1], "ζ")
    data = n.mdl.mdl.args
    y = n.chn[:year][1]
    return plotdensrbf(coefs, data, "Year $y", clim)
end

function plotdensrbf(coefs::Vector{Float64}, data::NamedTuple,
                     title = nothing, clim = nothing)
    cx = data.cx
    cy = data.cy
    scale = data.rbf_scale
    R = data.R
    return plotdensrbf_(coefs, cx, cy, scale, R, title, clim)
end

function plotdensrbf_(coefs::Vector{Float64},
                      cx::Vector{Float64},
                      cy::Vector{Float64},
                      scale::Float64,
                      R::Vector{Float64}, title = nothing, clim = nothing)
    Rmin, Rmax = extrema(R)
    vals = range(Rmin, Rmax, 1000)
    ## reversing y-values so we don't have to use heatmap(.., yflip =
    ## true), because this also flips the axis labels
    mat = [interpolant(rbf, Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale)
           for Rfrom in vals, Rto in vals]';
    mat = mat .- mean(mat)
    p = Plots.heatmap(vals, vals, mat, title = title, clim = clim, aspect_ratio = :equal)
    return mat, p
end
