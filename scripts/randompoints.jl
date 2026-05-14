using Distances, CairoMakie, Random, StatsBase, Distributions
include("/home/konstantin/code/src/plotutils.jl") ## helper functions for Makie
include("/home/konstantin/paper/plotting/utils.jl")

function generate_points(N)
    a, b = genab()
    x = rand(N) .* a
    y = rand(N) .* b
    return hcat(x, y)
end

function sequential(points, p = .1; start_idx = nothing)
    N = size(points, 1)
    i = isnothing(start_idx) ? rand(1:N) : start_idx

    dists = pairwise(Euclidean(), points', points[i:i, :]', dims=2)[:]
    order = sortperm(dists)
    order = filter(j -> j != i, order) # exclude starting point

    # sequential search
    for (k, j) in enumerate(order)
        if rand() < p
            return (dists[j], j, k)  # success
        end
    end
    return (NaN, nothing, length(order))
end

function radial(points, m, ϕ, i = nothing)
    i = isnothing(i) ? rand(1:size(points, 1)) : i
    ds = pairwise(Euclidean(), points', points[i:i, :]', dims=2)[:]
    ds = ds[ds .> 0]

    if rand() < ϕ
        return rand(ds)
    else
        ds = ds[ds .< m]
        if length(ds) > 0
            return rand(ds)
        else
            return NaN
        end
    end
end

function target_decline(points, γ = 2, start_idx = nothing)
    N = size(points, 1)
    i = isnothing(start_idx) ? rand(1:N) : start_idx
    dists = pairwise(Euclidean(), points', points[i:i, :]', dims=2)
    dists = dists[dists .> 0]
    return StatsBase.sample(dists, Weights(1 ./ dists .^γ))
end

function simulate(N, dist, P, p1 = .2, p2 = .5, ϕ = .1, γ = 2)
    points = generate_points(N)
    D = pairwise(Euclidean(), points')
    S = 10^3
    td = [target_decline(points, γ)[1] for _ in 1:S];
    sd = [sequential(points, p1)[1] for _ in 1:S];
    md = [sequential(generate_points(P), p2)[1] for _ in 1:S];
    rd = [radial(points, rand(dist), ϕ)[1] for _ in 1:S]

    lbls = ["Insensitive", "Target, γ = $γ", "Sequential, p = $p1",
            "Sequential, N = $P, p = $p2", "Radial, ϕ = $ϕ"]
    series = [vec(D[D .> 0]), td, sd, md[.!isnan.(md)], rd[.!isnan.(rd)]]
    return visualize(series, lbls)
end

function visualize(series, lbls, main, fig = genfig((10, 6)))
    ax = Axis(fig[1, 1], xlabel = "Distance (km)",
              ylabel = "CDF", title = main)
    xs = 1:821
    for (i, (s, l)) in enumerate(zip(series, lbls))
        lines!(ax, xs, ecdf(s).(xs), label = l, color = Cycled(i))
    end
    axislegend(ax; position = :rb)
    return fig
end

function simulate_radial(N, D, main = "", fn = nothing)
    points = generate_points(N)
    S = 10^3
    ϕ = [.0, .1, .15]
    td = [target_decline(generate_points(400), 2)[1] for _ in 1:S];
    series = [td]
    lbls = ["Target"]
    for p in ϕ
        rd = [radial(points, rand(D), p) for _ in 1:S]
        lbl = "ϕ = $p"
        push!(series, rd[.!isnan.(rd)])
        push!(lbls, lbl)
    end
    if !isnothing(fn); save(fn, fig); end
    return visualize(series, lbls, main)
end

function simulate_sequential(ns, ps, main = "", fn = nothing)
    S = 10^3
    td = [target_decline(generate_points(400), 2) for _ in 1:S]
    series = [td]
    lbls = ["Target"]
    for n in ns
        for p in ps
            points = generate_points(n)
            sam = [sequential(points, p)[1] for _ in 1:S]
            lbl = "N = $n, p = $p"
            push!(series, sam[.!isnan.(sam)])
            push!(lbls, lbl)
        end
    end
    if !isnothing(fn); save(fn, fig); end
    return visualize(series, lbls, main)
end

function genab()
    A = 357000
    a = sqrt(A / 1.33)
    b = 1.33a
    return a, b
end

function plotgamma(d, fn = nothing)
    fig = genfig((6, 5));
    ax = Axis(fig[1, 1], xlabel = "Distance (km)")
    lines!(ax, d)
    hideydecorations!(ax)
    if !isnothing(fn); save(fn, fig); end
    return fig
end

plotgamma(Gamma(2, 35/1))
simulate_radial(400, Gamma(2, 35/1), "Radial Search, Gamma(2, 35/1)")
simulate_sequential([10, 30], [.7], "Modified Sequential Search")
simulate_sequential([400], [.05, .1, .2], "Sequential Search")
