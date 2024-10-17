using Pkg
Pkg.activate(".")
using Turing, CSV, DataFrames, StatsPlots, Plots, LaTeXStrings

fit_model = function(df)
    @model function densmod(fd, td, y)
        a ~ Normal(0, 5)
        b ~ Normal(1, 5)
        s ~ truncated(Normal(3, 3), 0, Inf)

        preds = [exp(td[i] - fd[i]) ^ (a + b * exp(fd[i])) for i in 1:length(fd)]
        y ~ arraydist([Normal(m, s) for m in preds])
        return preds
    end

    model = densmod(df.fromdens, df.todens, df.funval)                       
    chain = Turing.sample(model, NUTS(), 1000)
    df.preds = generated_quantities(model, chain[1000])[1]
    plot(chain)
    return df
end

plot_surface = function(df, max)
    Plots.surface(df.todens, df.fromdens, df.funval,
                  title = "Cheby density function evaluation\n
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

df = CSV.read("./manuscript_input/germdensfun_18-25.csv", DataFrame)
df = df[df.fromdens .< 2.7, :]
df = df[df.todens .< 2.7, :]

df = fit_model(df)
plot_prediction(df)

histogram(df.funval)
histogram(df.fromdens)
scatter(df.fromdens, df.funval)
scatter(df.todens, df.funval)
sort!(df)


plot_surface(df, 2.7)



df[df.fromdens .< 2.7, :]
