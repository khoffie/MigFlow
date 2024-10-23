using Pkg
Pkg.activate(".")
using Turing, CSV, DataFrames, StatsPlots, Plots, LaTeXStrings, Serialization

plot_surface = function()
    plot_surface_ = function(df, age, max)
        Plots.surface(df.todens, df.fromdens, df.funval,
                      title = "Cheby density, agegroup $(age)\n
                  max logreldens = $(max)",
                      xlab = L"$\rho_o$", ylab = L"$\rho_d$",
                      camera = (10, 45))
    end
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    for age in ages
        df = CSV.read("./manuscript_input/germdensfun_$(age).csv", DataFrame)
        df = df[df.fromdens .< 2.7, :]
        df = df[df.todens .< 2.7, :]
        p = plot_surface_(df, age, 2.7)
        savefig("./manuscript_input/plots/dens_$(age).pdf")
    end
end

savelp = function()
    path = "./manuscript_input"
    out = Dict{String, Any}()  # Create a dictionary to store the plots
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    for age in ages
        chain = deserialize(path * "/germchain_$(age)")
        out[age] = Plots.plot(chain[:lp], title = "LP for $(age)")
    end
    p = Plots.plot(out["below18"], out["18-25"], out["25-30"],
                   out["30-50"], out["50-65"], out["above65"],
                   layout = (3, 2))
    savefig(path * "/plots/lp_values.pdf")
end

savelp()
plot_surface()
