
using Pkg
Pkg.activate("..")
using Turing, CSV, DataFrames, StatsPlots, Plots, LaTeXStrings
import GeoDataFrames as GDF

fit_model = function(df)
    @model function densmod(fd, td, y)
        a ~ Normal(0, 5)
        b ~ Normal(1, 5)
        c ~ Normal(1, 5)
        d ~ Normal(1, 5)
        s ~ truncated(Normal(3, 3), 0, Inf)

        preds = [a + b * td[i] + c * td[i]^2 + d * td[i]^3  for i in eachindex(fd)]
        y ~ arraydist([Normal(m, s) for m in preds])
        return preds
    end

    model = densmod(df.fromdens, df.todens, df.funval)                       
    chain = Turing.sample(model, NUTS(), 1000)
    df.preds = generated_quantities(model, chain[1000])[1]
    display(plot(chain))
    return df
end

fit_model = function(df)
    @model function densmod(fd, td, y)
        α_0 ~ Normal(0, 2)
        α_1 ~ Normal(1, 2)
        β_0 ~ Normal(0, 2)
        β_1 ~ Normal(1, 2)
        γ_0 ~ Normal(0, 2)
        γ_1 ~ Normal(1, 2)
        δ_0 ~ Normal(0, 2)
        δ_1 ~ Normal(1, 2)
        s ~ truncated(Normal(3, 3), 0, Inf)

        a = α_0 .+ α_1 .* fd
        b = β_0 .+ β_1 .* fd
        c = γ_0 .+ γ_1 .* fd
        d = δ_0 .+ δ_1 .* fd
##        preds = [a[i] + b[i] * td[i] + c[i] * td[i]^2 + d[i] * td[i]^3  for i in eachindex(fd)]
        preds = [a[i] + b[i] * td[i] + c[i] * td[i]^2  for i in eachindex(fd)]

        y ~ arraydist([Normal(m, s) for m in preds])
        return preds
    end

    model = densmod(df.fromdens, df.todens, df.funval)                       
    chain = Turing.sample(model, NUTS(), 1000)
    df.preds = generated_quantities(model, chain[1000])[1]
    display(plot(chain))
    return df
end

plot_surface = function(df, age, max)
    Plots.surface(df.todens, df.fromdens, df.funval,
                  title = "Cheby density, agegroup $(age)\n
                  max logreldens = $(max)",
                  xlab = L"$\rho_o$", ylab = L"$\rho_d$",
                  camera = (10, 45))
end

plot_prediction = function(df)
    cc = round(cor(df.preds, df.funval), digits = 2)
    Plots.plot(df.preds, df.funval, seriestype = :scatter,
         xlabel = "Prediction", ylabel = "Cheby value",
         title = "Density model", label = "Correlation = $(cc)")
end

ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
for age in ages
    df = CSV.read("../manuscript_input/germdensfun_$(age).csv", DataFrame)
    df = df[df.fromdens .< 2.7, :]
    df = df[df.todens .< 2.7, :]
    p = plot_surface(df, age, 2.7)
	display(p)
end

df = CSV.read("../manuscript_input/germdensfun_18-25.csv", DataFrame)
df = df[df.fromdens.<2.7, :]
df = df[df.todens.<2.7, :]
districts = CSV.read("../data/districts.csv", DataFrame)
md = median(districts.density)
districts.logreldens = log.(districts.density ./ md)
f_shp = "/home/konstantin/Diss/inst/extdata/clean/shapes/districts.shp"
shp = GDF.read(f_shp)
shp[shp.geometry != "Geometry: wkbPolygon", :]
filter(:geometry => ==("Geometry: wkbPolygon"), shp)
buffer.(shp.geometry, 10)
plot(shp.geometry[1:2], fill_z = )

shp.AGS = parse.(Int64, shp.AGS)
to_join = select(districts, [:distcode, :logreldens])
shp = innerjoin(shp,  to_join, on = :AGS => :distcode)

test = shp.geometry[3]
push!(test

shp.geometry
println(length(shp.geometry))  # Number of geometries (districts)
println(length(shp.logreldens))  # Number of logreldens values

plot(shp.geometry, fill_z = shp.logreldens)
plot(shp.geometry[3])
df = fit_model(df)
plot_prediction(df)

plot_surface(df, "18-25", 2.7)

histogram(df.funval)
histogram(df.fromdens)
histogram(df.todens, xlab = "todens")

p1 = scatter(df.fromdens, df.funval, xlab = "fromdens", ylab = "cheby val")
p2 = scatter(df.todens, df.funval, xlab = "todens", ylab = "cheby val")

p1 = histogram(districts.density, title = "Histogram Density")
p2 = histogram(districts.logreldens, title = "Histogram Logreldens")
plot(p1, p2)

