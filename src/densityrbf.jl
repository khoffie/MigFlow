function plotdtf(r::EstimationResult, crange = nothing, fig = genfig((10, 10)),
                 x = 1, y = 1, legend = true)
    m, a, yr = dtfmat(r)
    if isnothing(crange); crange = calc_crange(m); end
    plotdtf(m, a, yr, crange, fig, x, y, legend)
    return m, fig
end

function plotdtf(m::LinearAlgebra.Adjoint{Float64, Matrix{Float64}},
                 a::String, yr::Int, crange, fig, x, y, legend = true)
    qs = [.01, .25, .5, .75, .9, 1]
    ## from round.(quantile(districts.density, [0.0, .25, .5, .75, 0.9, 1.0]), digits = 0) in 2017
    ys = string.([35, 115, 200, 670, 1500, 4700])
    tks = (100qs, ys)

    ax = Axis(fig[x, y],
              xlabel = L"\textrm{Origin} $\mathrm{pop} / \mathrm{km}^2$",
              ylabel = L"\textrm{Destination} $\mathrm{pop} / \mathrm{km}^2$",
              xticks = tks, yticks = tks,
              title = string(yr), aspect = DataAspect(),
              xticklabelrotation = 1.0)
    hm = Makie.heatmap!(ax, m, colorrange = crange, colormap = :roma)
    if legend
        Colorbar(fig[:, y + 1], hm)
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
    anchor = data.anchor
    Rmin = -1.0; Rmax = 1.0
    vals = range(Rmin, Rmax, 100)
    m = [interp_anchor(Rfrom, Rto, coefmat(coefs ./ 10), cx, cy, scale, anchor)
         for Rto in vals, Rfrom in vals]'
    return m, a, y
end
