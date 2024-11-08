plot_surface = function(path)
    plot_surface_ = function(df, age, max)
        Plots.surface(df.fromdens, df.todens, df.funval,
                      title = "Cheby density, agegroup $(age)\n
                  max logreldens = $(max)",
                      xlab = L"$\rho_o$", ylab = L"$\rho_d$",
                      camera = (45, 45))
    end
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    years = 2011 : 2017
    for y in years
        for a in ages
            df = CSV.read(joinpath(path,  "germdensfun_$(y)_$(a).csv"), DataFrame)
            df = df[df.fromdens .< 2.7, :]
            df = df[df.todens .< 2.7, :]
            p = plot_surface_(df, a, 2.7)
            savefig(joinpath(path, "plots/dens_$(y)_$(a).pdf"))
        end
    end
end

savelp = function(path, from = nothing, to = nothing)
    function plotlp(chain, age, from = nothing, to = nothing)
        ss = length(chain)
        # lastlp = chain[:lp][ss, :].data
        # lastlp = round(mean(lastlp), digits = 0)
        from = isnothing(from) ? 1 : from
        to = isnothing(to) ? ss : to
        xvals = [from:to]
        p = Plots.plot(xvals, chain[:lp][from : to, :],
                       title = "LP for $(age)", label = "")
        return p
    end
    out = Dict{String, Any}()  # Create a dictionary to store the plots
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    years = 2011 : 2017
    for y in years
        for a in ages
            chain = deserialize(joinpath(path ,"germchain_$(y)_$(a).csv"))
            out[a] = plotlp(chain, a, from, to)
        end
        p = Plots.plot(out["below18"], out["18-25"], out["25-30"],
                       out["30-50"], out["50-65"], out["above65"],
                       layout = (3, 2))
    display(p)
    savefig(joinpath(path, "plots/lp_values_$(y).pdf"))
    end
end

plot_distance = function(path)
    plot_distance_ = function(df, age)
        Plots.plot(df.dist, df.out, seriestype = :scatter, alpha = .1,
                   xlab = "Distance", ylab = "log.((df.flows .+ 1) ./ df.preds)",
                   title = "agegroup $(age)")
    end
    out = Dict{String, Any}()  # Create a dictionary to store the plots
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    years = 2011 : 2017
    for y in years
        for a in ages
            df = CSV.read(joinpath(path, "germflows_$(y)_$(a).csv"), DataFrame)
            N = nrow(df)
            size = 2000
            df = df[Random.shuffle(1:N)[1:min(size, N)], :]
            df.out = log.((df.flows .+ 1) ./ df.preds)
            out[a] = plot_distance_(df, a)
        end
        p = Plots.plot(out["below18"], out["18-25"], out["25-30"],
                       out["30-50"], out["50-65"], out["above65"],
                       layout = (3, 2))
        display(p)
        savefig(joinpath(path,  "plots/distance_$(y).pdf"))
    end
end

post_process = function(path = nothing, lp_from = nothing, lp_to = nothing, render_doc = true)
    isnothing(path) ? path = path = readline("./writeup/juliaout_path.txt") : path
    println(path)
    mkpath(joinpath(path, "plots"))
    savelp(path, lp_from, lp_to)
    plot_surface(path)
    plot_distance(path)
    println("Plots saved")

    # report.Rmd reads julia_output_path from file = "./writeup/juliaout_path.txt"
    f = "./writeup/_main.Rmd"
    if isfile(f);  rm(f); end
    render_doc ? R"helpeR::render_doc('./writeup', 'report.Rmd')" : println("Report not generated")
end
