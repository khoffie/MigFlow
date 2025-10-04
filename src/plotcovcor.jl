function plotcovcor(r::EstimationResult)
    a, y = getageyear(r)
    crm =  cov2cor(r.cvm)
    nms = r.ses.name
    n = length(nms)
    fig = Figure(size = (1600, 800));
    ax1 = Axis(fig[1, 1],
              title = "Correlation Matrix",
              xticks = (1:n, nms),
              yticks = (1:n, nms)
              #xticklabelrotation = π/2,  # rotate to avoid overlap
              )
    hm = heatmap!(ax1, Matrix(crm))
    Colorbar(fig[2, 1], hm, label = "Correlation", vertical = false)

    ax2 = Axis(fig[1, 2],
              title = "Covariance Matrix",
              xticks = (1:n, nms),
              yticks = (1:n, nms)
              #xticklabelrotation = π/2,  # rotate to avoid overlap
              )
    hm = heatmap!(ax2, Matrix(r.cvm))
    Colorbar(fig[2, 2], hm, label = "Covariance", vertical = false)
    prettytitle!(fig, "$a, $y")
    return fig
end
