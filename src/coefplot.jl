function coefplot(r::EstimationResult)
    a, y = getageyear(r)
    fig = Figure(size = (1000, 400), fontsize = 10);
    ax = Axis(fig[1, 1],
              xlabel = "Parameter",
              ylabel = "Estimate +. 2se",
              xticks = (1:nrow(ses), ses.name),
              xgridvisible = false)
    errorbars!(ax, 1:nrow(ses), ses.coef, 2ses.se)
    Makie.scatter!(ax, 1:nrow(ses), ses.coef)
    prettytitle!(fig, "Estimates $a in $y")
    return fig
end
