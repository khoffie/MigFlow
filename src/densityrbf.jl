function plotdtf(r::EstimationResult, crange = nothing,
                 fig = Figure(size=(400, 400), fontsize = 10),
                 x = 1, y = 1, legend = true)
    m, a, yr = dtfmat(r)
    if isnothing(crange); crange = extrema(m); end
    plotdtf(m, a, yr, crange, fig, x, y, legend)
    return m, fig
end

function plotdtf(m::Matrix{Float64}, age::String, year::Int, crange,
                 fig, x, y, legend = true)
    ax = Axis(fig[x, y],
              xlabel = "Origin density",
              ylabel = "Destination density",
              xticks = ([1, 100], ["low", "high"]),
              yticks = ([1, 100], ["low", "high"]),
              title = string(year), aspect = DataAspect())
    hm = Makie.heatmap!(ax, m, colorrange = crange)
    if legend
        Colorbar(fig[:, y + 1], hm, vertical = true, width = 5)
        prettytitle!(fig, "Density Transition Function, $age")
    end
end

function dtfmat(r::EstimationResult)
    a, y = getageyear(r)
    coefs = extract_coefs(r.chn[end, :, 1], "ζ")
    data = r.mdl.mdl.args
    cx = data.cx
    cy = data.cy
    scale = data.rbf_scale
    R = data.R
    Rmin, Rmax = extrema(R)
    anchor = median(R)
    vals = range(Rmin, Rmax, 100)
    m = [interp_anchor(Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale, anchor)
         for Rto in vals, Rfrom in vals]';
    return m .- mean(m), a, y
end
