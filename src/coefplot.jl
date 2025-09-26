function coefplot(r::EstimationResult)
    df = r.ses
    xs = 1:nrow(df)
    a, y = getageyear(r)
    fig = Figure(size = (1000, 400), fontsize = 10);
    ax = Axis(fig[1, 1],
              xlabel = "Parameter",
              ylabel = "Estimate +- 2se",
              xticks = (xs, df.name),
              xgridvisible = false)
    errorbars!(ax, xs, df.coef, 2df.se)
    Makie.scatter!(ax, xs, df.coef)
    prettytitle!(fig, "Estimates $a in $y")
    return fig
end
