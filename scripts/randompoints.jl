using Distances, CairoMakie, Random, StatsBase, Distributions
include("/home/konstantin/code/src/plotutils.jl") ## helper functions for Makie
include("/home/konstantin/paper/plotting/utils.jl")

function generate_points(N, a, b)
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
    points = generate_points(N, a, b)
    D = pairwise(Euclidean(), points')
    S = 10^3
    sd = [sequential(points, p1)[1] for _ in 1:S];
    md = [sequential(generate_points(P, a, b), p2)[1] for _ in 1:S];
    md = md[.!isnan.(md)]
    td = [target_decline(points, γ)[1] for _ in 1:S];


    fig = genfig();
    ax = Axis(fig[1, 1], xlabel = "Distance (km)",
              ylabel = "CDF")
    xs = 1:821
    lbls = ["Insensitive", "Target, γ = $γ", "Sequential, p = $p1",
            "Sequential, N = $P, p = $p2", "Radial, ϕ = $ϕ"]
    series = [vec(D[D .> 0]), td, sd, md, rd[.!isnan.(rd)]]
    for (i, (s, l)) in enumerate(zip(series, lbls))
        lines!(ax, xs, ecdf(s).(xs), label = l, color = Cycled(i))
    end
    axislegend(ax; position = :rb)

    fig2 = genfig();
    ax = Axis(fig2[1, 1], aspect = DataAspect(), title = "Germany")
    scatter!(ax, generate_points(80, a, b))

    return (; d = fig, g = fig2)
end

N = 10^3
A = 357000
a = sqrt(A / 1.33)
b = 1.33a
lines(Gamma(5, 40/4))
simulate(400, Gamma(5, 40/4), 10, .2, .9, .15, 2).d
