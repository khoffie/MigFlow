function plotdtf(r::EstimationResult, crange = nothing,
                 fig = Figure(size=(400, 400), fontsize = 10),
                 x = 1, y = 1, legend = true)
    m, a, yr = dtfmat(r)
    if isnothing(crange); crange = extrema(m); end
    plotdtf(m, a, yr, crange, fig, x, y, nothing, legend)
    return m, fig
end

function plotdtf(r::EstimationResult, districts::DataFrame, crange = nothing,
                 fig = Figure(size=(400, 400), fontsize = 10),
                 x = 1, y = 1, legend = true)
    m, a, yr = dtfmat(r)
    if isnothing(crange); crange = extrema(m); end
    plotdtf(m, a, yr, crange, fig, x, y, districts, legend)
    return m, fig
end

function plotdtf(m::LinearAlgebra.Adjoint{Float64, Matrix{Float64}},
                 a::String, yr::Int, crange, fig, x, y,
                 districts = nothing, legend = true)
    if !isnothing(districts)
        di = year(districts, yr)
        l, h = Int.(round.(extrema(di.density), digits = 0))
        l = string(l); h = string(h)
        xt = ([1, 100], [l, h])
        yt = ([1, 100], [l, h])
    else
        xt = ([1, 100], ["low", "high"])
        yt = ([1, 100], ["low", "high"])
    end
    ax = Axis(fig[x, y],
              xlabel = "Origin density",
              ylabel = "Destination density",
              xticks = xt,
              yticks = yt,
              title = string(yr), aspect = DataAspect())
    hm = Makie.heatmap!(ax, m, colorrange = crange)
    if legend
        Colorbar(fig[:, y + 1], hm, vertical = true, width = 5)
        prettytitle!(fig, "Density Transition Function, $a")
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
    Rmin, Rmax = data.Rmin, data.Rmax
    anchor = data.anchor
    vals = range(Rmin, Rmax, 100)
    m = [interp_anchor(Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale, anchor)
         for Rto in vals, Rfrom in vals]'
    return m, a, y
end
