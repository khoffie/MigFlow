using Pkg
Pkg.activate(".")
using Turing, CSV, DataFrames, StatsPlots, Plots, LaTeXStrings, Serialization, Random

plot_surface = function(path)
    plot_surface_ = function(df, age, max)
        Plots.surface(df.fromdens, df.todens, df.funval,
                      title = "Cheby density, agegroup $(age)\n
                  max logreldens = $(max)",
                      xlab = L"$\rho_o$", ylab = L"$\rho_d$",
                      camera = (45, 45))
    end
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    for age in ages
        df = CSV.read(joinpath(path,  "germdensfun_$(age).csv"), DataFrame)
        df = df[df.fromdens .< 2.7, :]
        df = df[df.todens .< 2.7, :]
        p = plot_surface_(df, age, 2.7)
        savefig(joinpath(path, "plots/dens_$(age).pdf"))
    end
end

savelp = function(path, from = nothing, to = nothing)
    out = Dict{String, Any}()  # Create a dictionary to store the plots
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    for age in ages
        chain = deserialize(joinpath(path ,"germchain_$(age).csv"))
        ss = length(chain)
        lastlp = chain[:lp][ss]
        lastlp = round(lastlp, digits = 0)
        from = isnothing(from) ? 1 : from
        to = isnothing(to) ? ss : to
        xvals = [from:to]
        out[age] = Plots.plot(xvals, chain[:lp][from:to],
                              title = "LP for $(age)",
                              label = "lastlp: $(lastlp)",
                              xlab = "Sample", ylab = "LP")
    end
    p = Plots.plot(out["below18"], out["18-25"], out["25-30"],
                   out["30-50"], out["50-65"], out["above65"],
                   layout = (3, 2))
    savefig(joinpath(path, "plots/lp_values.pdf"))
end

function plotlp(chain, age, from = nothing, to = nothing)
    ss = length(chain)
    # lastlp = chain[:lp][ss, :].data
    # lastlp = round(mean(lastlp), digits = 0)
    from = isnothing(from) ? 1 : from
    to = isnothing(to) ? ss : to
    xvals = [from:to]
    plt = Plots.plot(xvals, chain[:lp][from:to, :],
                     title = "LP for $(age)",
                     xlab = "Sample", ylab = "LP")
    return plt
end

plot_distance = function(path)
    plot_distance_ = function(df, age)
        Plots.plot(df.dist, df.out, seriestype = :scatter, alpha = .1,
                   xlab = "Distance", ylab = "log.((df.flows .+ 1) ./ df.preds)",
                   title = "agegroup $(age)")
    end
    out = Dict{String, Any}()  # Create a dictionary to store the plots
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    for age in ages
        df = CSV.read(joinpath(path, "germflows_$(age).csv"), DataFrame)
        N = nrow(df)
        size = 2000
        df = df[Random.shuffle(1:N)[1:min(size, N)], :]
        df.out = log.((df.flows .+ 1) ./ df.preds)
        out[age] = plot_distance_(df, age)
    end
    p = Plots.plot(out["below18"], out["18-25"], out["25-30"],
                   out["30-50"], out["50-65"], out["above65"],
                   layout = (3, 2))
    display(p)
    savefig(joinpath(path,  "plots/distance.pdf"))
end

post_process = function(path)
    mkpath(joinpath(path, "plots"))
    savelp(path)
    plot_surface(path)
    Plots.scalefontsizes(.4)    
    plot_distance(path)
    println("Plots saved")

    file = "./writeup/juliaout_path.txt"
    open(file, "w") do f
        write(f, path)
    end
    # report.Rmd reads julia_output_path from file = "./writeup/juliaout_path.txt"
    R"helpeR::render_doc('~/Documents/GermanMigration/writeup', 'report.Rmd')"
    R"helpeR::render_doc('~/Documents/GermanMigration/writeup', 'math.tex')"
    println("Report generated")
end


chain = deserialize(joinpath("manuscript_input/2024-10-29_16-20-27", "germchain_25-30.csv"))
