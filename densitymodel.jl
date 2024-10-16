using Turing, CSV, DataFrames, StatsPlots

@model function densmod(fd, td, y)
    a ~ Normal(0, 5)
    b ~ Normal(1, 5)
    s ~ truncated(Normal(3, 3), 0, Inf)

    for i in 1:length(y)
        ratio = exp(td[i] - fd[i])
        m = ratio ^ (a + b * exp(fd[i]))
        y[i] ~ Normal(m , s )
    end
end

gen_preds = function(df, chain)
    predict = function(fd, td, a, b)
        p = exp(td - fd) ^ (a + b * exp(fd))
        return p
    end

    df.preds = zeros(nrow(df))
    for i in 1:nrow(df)
        df.preds[i] = predict(df.fromdens[i], df.todens[i], mean(chain[:a]), mean(chain[:b]))
    end
    return df
end

df = CSV.read("./manuscript_input/germdensfun_18-25.csv", DataFrame)
histogram(df.funval)
scatter(df.fromdens, df.funval)
scatter(df.todens, df.funval)


model = densmod(df.fromdens, df.todens, df.funval)                       
chain = Turing.sample(model, NUTS(), 1000)


gen_preds(df, chain)
    

scatter(df.preds, df.funval)

