using Distances, CairoMakie, Random, StatsBase, Distributions
include("/home/konstantin/code/src/plotutils.jl") ## helper functions for Makie
include("/home/konstantin/paper/plotting/utils.jl")

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

function simulate(N, p = .2)
    points = generate_points(N, a, b)
    D = pairwise(Euclidean(), points')

    md = [simulate_move(points, p)[1] for _ in 1:10^3];
    md = md[.!isnan.(md)]
    qd = [quadratic_decline(points)[1] for _ in 1:10^3];

    fig = genfig();
    ax = Axis(fig[1, 1], xlabel = "Distance (km)",
              title = "Number points = $N, success prob = $p",
              yticks = cdfticks(), ylabel = "CDF")
    xs = 1:821
    lbls = ["Insensitive", "Quadratic", "Sequential"]
    series = [vec(D[D .> 0]), qd, md]
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
res = simulate(N)
res.d

function simulate_radial(points, m, i = nothing)
    i = isnothing(i) ? rand(1:size(points, 1)) : i
    ds = pairwise(Euclidean(), points', points[i:i, :]', dims=2)[:]
    ds = ds[ds .> 0 .&& ds .< m]
    if length(ds) > 0
        return rand(ds)
    else
        return NaN
    end
end

m = 7000


moves = zeros(10^3)
for i in 1:10^3
    m = rand(Gamma(5, 30/4))
    moves[i] = simulate_radial(points, m)
end

e = ecdf(moves[.!isnan.(moves)])
