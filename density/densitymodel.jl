using Pkg
Pkg.activate("..")
using Turing, CSV, DataFrames, StatsPlots, Plots, LaTeXStrings
import GeoDataFrames as GDF

plot_surface = function(df, age, max)
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
    p = plot_surface(df, age, 2.7)
    savefig("./manuscript_input/plots/dens_$(age).pdf")
end
