function dtfmat(r::EstimationResult)
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

function plotdtf(r::EstimationResult)
    m, a, y = dtfmat(r)
    fig = Figure(size=(400, 400), fontsize = 10)
    plotdtf(m, a, y, extrema(m), fig, 1, 1)
    return fig
end

function plotdtf(m::Matrix{Float64}, age::String, year::Int, crange, f, x, y)
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
