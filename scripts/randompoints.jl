using Distances, CairoMakie, Random, StatsBase
include("/home/konstantin/code/src/plotutils.jl") ## helper functions for Makie

function generate_points(N, a, b)
    x = rand(N) .* a
    y = rand(N) .* b
    return hcat(x, y)
end

function simulate_move(points, p = .1; start_idx = nothing)
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

function quadratic_decline(points; start_idx = nothing)
    N = size(points, 1)
    i = isnothing(start_idx) ? rand(1:N) : start_idx
    dists = pairwise(Euclidean(), points', points[i:i, :]', dims=2)
    dists = dists[dists .> 0]
    return StatsBase.sample(dists, Weights(1 ./ dists .^2))
end


N = 100

A = 357000
a = sqrt(A / 1.33)
b = 1.33a



function simulate(N, p = .2)
    points = generate_points(N, a, b)
    D = pairwise(Euclidean(), points', dims=2)

    md = [simulate_move(points, p)[1] for _ in 1:10^3];
    md = md[.!isnan.(md)]
    qd = [quadratic_decline(points)[1] for _ in 1:10^3];

    fig = genfig();
    ax = Axis(fig[1, 1], xlabel = "Distance (km)",
              title = "Number points = $N, success prob = $p")
    density!(ax, vec(D[D .> 0]), color = :blue, label = "All")
    density!(ax, qd, color = :green, label = "Quadratic")
    density!(ax, md, color = :red, label = "Sequential")

    axislegend(ax; position = :rt)

    fig2 = genfig();
    ax = Axis(fig2[1, 1], aspect = DataAspect(), title = "Germany")
    scatter!(ax, generate_points(80, a, b))

    return (; d = fig, g = fig2)
end


simulate(150, .07)[1]
