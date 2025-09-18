# function plotdensrbf(n::NamedTuple, clim)
#     coefs = extract_coefs(n.chn[end, :, 1], "ζ")
#     data = n.mdl.mdl.args
#     y = n.chn[:year][1]
#     return plotdensrbf(coefs, data, "Year $y", clim)
# end

# function plotdensrbf(coefs::Vector{Float64}, data::NamedTuple,
#                      title = nothing, clim = nothing)
#     cx = data.cx
#     cy = data.cy
#     scale = data.rbf_scale
#     R = data.R
#     return plotdensrbf_(coefs, cx, cy, scale, R, title, clim)
# end

# function plotdensrbf_(coefs::Vector{Float64},
#                       cx::Vector{Float64},
#                       cy::Vector{Float64},
#                       scale::Float64,
#                       R::Vector{Float64}, title = nothing, clim = nothing)
#     Rmin, Rmax = extrema(R)
#     vals = range(Rmin, Rmax, 1000)
#     ## reversing y-values so we don't have to use heatmap(.., yflip =
#     ## true), because this also flips the axis labels
#     mat = [interpolant(rbf, Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale)
#            for Rfrom in vals, Rto in vals]';
#     mat = mat .- mean(mat)
#     p = Plots.heatmap(vals, vals, mat, title = title, clim = clim, aspect_ratio = :equal)
#     return mat, p
# end

function dtfmat(r)
    a, y = getageyear(r)
    coefs = extract_coefs(r.chn[end, :, 1], "ζ")
    data = r.mdl.mdl.args
    cx = data.cx
    cy = data.cy
    scale = data.rbf_scale
    R = data.R
    Rmin, Rmax = extrema(R)
    vals = range(Rmin, Rmax, 100)
    m = [interpolant(rbf, Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale)
         for Rto in vals, Rfrom in vals]';
    return m .- mean(m), a, y
end

function plotdtf(m, age, year, crange, f, x, y)
    ax = Axis(f[x, y],
              xlabel = "Origin density",
              ylabel = "Destination density",
              xticks = ([1, 1000], ["low", "high"]),
              yticks = ([1, 1000], ["low", "high"]),
              title = string(year), aspect = DataAspect())
    Makie.heatmap!(ax, m, colorrange = crange)
end

function plotdtfyears(results, age, years, f, crange = nothing)
    o = loopstruct(results, dtfmat, age, years)
    if isnothing(crange)
        crange = extrema(reduce(vcat, [o[i][1] for i in eachindex(o)]))
    end
#     ncols = length(years) > 9 ? 3 : 2
    for i in eachindex(years)
        pos = grid_position(i, 3)
        m, a, y = dtfmat(results[age[1]][years[i]])
        plotdtf(m, a, y, crange, f, pos[1], pos[2])
    end
    Colorbar(f[:, end + 1], colorrange = crange)
    a = replace(string(age[1]), "age" => "")
    Label(f[0, :], "Density transition function, $a",
          fontsize = 20, tellwidth = false)
    return f
end
