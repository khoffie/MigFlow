using Pkg
Pkg.activate(".")
using Turing, CSV, DataFrames, StatsPlots, Plots

df = CSV.read("./manuscript_input/germdensfun_18-25.csv", DataFrame)
# districts = CSV.read("./data/districts.csv", DataFrame)
histogram(df.funval)
scatter(df.fromdens, df.funval)
scatter(df.todens, df.funval)

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
plot(chain)
df.preds = generated_quantities(model, chain[1000])[1]

cc = round(cor(df.preds, df.funval), digits = 2)
plot(df.preds, df.funval, seriestype = :scatter,
     xlabel = "Prediction", ylabel = "Cheby value",
     title = "Density model", label = "Correlation = $(cc)")


