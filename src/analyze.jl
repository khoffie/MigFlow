abstract type AbstractDataSet end

struct Flows <: AbstractDataSet
    df::DataFrame
end

struct Net <: AbstractDataSet
    df::DataFrame
end

struct AnalysisResult
    df::Flows ## flows df with predictions
    net::Net ## net, nmr, asymmetries
    quick::DataFrame ## deviance, MAE, MAE0, Skillscore
    asym::DataFrame ## bivariate asymmetries
    fig::Figure ## Main analysis plot
end

function Base.vcat(x::T, y::T) where {T<:AbstractDataSet}
    T(vcat(x.df, y.df))
end

function DataFrames.filter(f, x::T; kwargs...) where {T<:AbstractDataSet}
    T(DataFrames.filter(f, x.df; kwargs...))
end

DataFrames.groupby(x::AbstractDataSet, cols) = DataFrames.groupby(x.df, cols)
Base.names(x::AbstractDataSet) = Base.names(x.df)
Base.sort(x::AbstractDataSet, c) = Base.sort(x.df, c)

function analyze(r::EstimationResult, fig = genfig((20, 6)))
    df = modeldf(r)
    net = netdf(r)
    quick = quickdf(r)
    asym = asymdf(df)
    pointsize = 6
    ax1 = Axis(fig[1, 1],
               xlabel = L"\log \hat{y}",
               ylabel = L"\log y",
               title = L"\text{Mean deviance:}%$(quick.deviance[1])",
               aspect = DataAspect(),
               xgridvisible = false, ygridvisible = false)
    plotfit!(ax1, df.df.flows, df.df.preds, pointsize)

    tks = ([-1.0, -.5, 0.0, .5, 1.0], ["-1", "-.5", "0", ".5", "1"])
    ax2 = Axis(fig[1, 2],
               xlabel = L"(\hat{i} - \hat{o}) / (\hat{i} + \hat{o})",
               ylabel = L"(i - o) / (i + o)",
               title = L"\textrm{skillscore} = %$(round(quick.skillscore[1], digits = 2))",
               aspect = DataAspect(),
               xgridvisible = false, ygridvisible = false, xticks = tks, yticks = tks)
    Makie.ylims!(ax2, (-1, 1))
    Makie.xlims!(ax2, (-1, 1))
    plotasym!(ax2, net, pointsize)

    tks = ([0, 200, 400, 600, 800], string.([0, 2, 4, 6, 8]))
    ax3 = Axis(fig[1, 3],
               xlabel = L"\text{Distance (100km)}",
               ylabel = L"\log(y / \hat{y})",
               xgridvisible = false, ygridvisible = false, xticks = tks)
##    ylims!(ax3, -2, 2)
    plotdist!(ax3, df.df.flows, df.df.preds, df.df.dist, pointsize)

    ax4 = Axis(fig[1, 4],
               xlabel = L"\log(A_o  P_d)",
               ylabel = L"\log(y / \hat{y})",
               xgridvisible = false, ygridvisible = false)
    plotpop!(ax4, df.df.flows, df.df.preds, df.df.A, df.df.P, pointsize)
    println(typeof(df))
    return AnalysisResult(df, net, quick, asym, fig)
end

function modeldf(r::EstimationResult)
    a, y = getageyear(r)
    data = r.mdl.mdl.args
    df = DataFrame(
        fromdist = data.from,
        todist = data.to,
        flows = data.Y,
        preds = r.prd,
        dist = 100data.D, ## scaling back to original, better grab ds?
        A = data.A,
        P = exp.(data.P[data.to]) ## bec log(P) is saved
    )
    return Flows(df)
end

function quickdf(r::EstimationResult)
    m, a, y = getmeta(r)
    df = modeldf(r).df
    net = netdf(r).df
    dev = round2(deviance2(df.flows, df.preds))
    err = 100mae(net.asyma, net.asymap)
    trivial = 100mae(net.asyma, 0)
    skill = skillscore(net.asyma, net.asymap)
    quick = DataFrame(model = m, agegroup = a, year = y, deviance = dev,
                      mae = err, mae0 = trivial, skillscore = skill)
    return quick
end


function netdf(r::EstimationResult)
    return addmeta(r, addnetcols(calcnet(r)))
end

function addmeta(r::EstimationResult, df::T) where {T<:AbstractDataSet}
    df = df.df
    m, a, y = getmeta(r)
    df.agegroup .= a
    df.year .= y
    df.model .= m
    first = ["model","agegroup", "year"]
    last = setdiff(names(df), first)
    return T(select(df, vcat(first, last)))
end

function calcnet(r::EstimationResult)
    df = modeldf(r)
    dfout = combine(DataFrames.groupby(df, [:fromdist]),
                    :flows => sum => :outflux,
                    :preds => sum => :outfluxp)
    dfin = combine(DataFrames.groupby(df, [:todist]),
                   :flows => sum => :influx,
                   :preds => sum => :influxp)
    net = innerjoin(dfout, dfin, on = :fromdist => :todist)
    net = innerjoin(net, unique(df.df, :fromdist)[!, [:fromdist, :A]], on = :fromdist)
    return Net(rename!(net, :fromdist => :lc))
end

function addnetcols(df::Net)
    net = df.df
    net.net = net.influx .- net.outflux
    net.total = net.influx .+ net.outflux
    net.asyma = net.net ./ net.total
    net.nmra = net.net ./ net.A
    net.netp = net.influxp .- net.outfluxp
    net.totalp = net.influxp .+ net.outfluxp
    net.asymap = net.netp ./ net.totalp
    net.nmrap = net.netp ./ net.A
    return Net(net)
end

subset(x, n) = StatsBase.sample(1:length(x), n)

function extract_coefs(chn::Chains)
    return vec(chn.value.data)
end

function extract_coefs(chn::Chains, string::String)
    nms = String.(names(chn))
    return chn.value[occursin.(string, nms)].data
end

function extract_sample(chn, type = "best")
    if type == "best"
        m = argmax(chn[:lp].data)
        chn = chn[m[1], :, m[2]]
    end
    return chn
end

function plotfit!(ax, flows, preds, size)
    function sub(df, N)
        idx = subset(flows, N)
        x = log.(preds)[idx]
        y = log.(flows)[idx]
        return DataFrame(; x, y)
    end
    df = sub(flows, 10^3)
    Makie.scatter!(ax, df.x, df.y, alpha = .5, markersize = size)
    diagonal!(ax, df.x, df.y)
    df = sub(flows, 10^4)
    smoother!(ax, df.x, df.y)
end

function plotasym!(ax, net::Net, size)
    net = net.df
    Makie.scatter!(ax, net.asymap, net.asyma, alpha = .5, markersize = size)
    diagonal!(ax, net.asymap, net.asyma)
    smoother!(ax, net.asymap, net.asyma)
end

function plotdist!(ax, flows, preds, dist, size)
    function sub(df, N)
        idx = subset(flows, N)
        y = pearres.(flows, preds)[idx]
##        y = log.(flows ./ preds)[idx]
        x = dist[idx]
        return sort(DataFrame(; x, y), :x)
    end
    df = sub(flows, 10^3)
    Makie.scatter!(ax, df.x, df.y, alpha = .5, markersize = size)
    Makie.hlines!(ax, [0], color = :darkred, linewidth = 2)
    df = sub(flows, 10^4)
    smoother!(ax, df.x, df.y)
end

function plotpop!(ax, flows, preds, frompop, topop, size)
    function sub(df, N)
        idx = subset(flows, 10^3)
        x = (log.(frompop) .+ log.(topop))[idx]
        y = (log.(flows ./ preds))[idx]
        return sort(DataFrame(; x, y), :x)
    end
    df = sub(flows, 10^3)
    Makie.scatter!(ax, df.x, df.y, alpha = .5, markersize = size)
    Makie.hlines!(ax, [0], color = :darkred, linewidth = 2)
    df = sub(flows, 10^4)
    smoother!(ax, df.x, df.y)
end

function deviance2(y, p)
    loss = zeros(length(y))
    for i in eachindex(y)
        loss[i] = y[i] * log(y[i] / p[i]) - (y[i] - p[i])
    end
    return 2mean(loss)
end

# mse(y, p) = mean((y .- p) .^ 2)
mae(y, p) = mean(abs.(y .- p))
## skillscore(y, p) = 1 - (mse(y, p) / mse(y, 0))
skillscore(y, p) = 1 - (mae(y, p) / mae(y, 0))
multires(y, p) = log(y / p)
devres(y, p) = sqrt((y * log(y / p)) - (y - p))
pearres(y, p) = (y - p) / sqrt(p)

function asymdf(df::Flows)
    df = df.df
    dfod = select(df, :fromdist, :todist, :flows => :outflux,
                  :preds => :outpreds)
    dfdo = select(df, :fromdist => :todist, :todist => :fromdist,
                  :flows => :influx, :preds => :inpreds)
    df1 = leftjoin(dfod, dfdo, on = [:fromdist, :todist])

    df1.total = df1.influx .+ df1.outflux
    df1.totalp = df1.inpreds .+ df1.outpreds
    df1.asyma = (df1.influx .- df1.outflux) ./ df1.total
    df1.asymap = (df1.inpreds .- df1.outpreds) ./ df1.totalp
    return dropmissing!(df1)
end

function getmeta(r::EstimationResult)
    a, y = getageyear(r)
    m = getmodel(r)
    return (; m, a, y)
end
