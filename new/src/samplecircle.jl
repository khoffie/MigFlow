using CSV, Plots, DataFrames, StatsBase, StatsPlots, Random
using Distributions, LazySets, Turing
include("src/utils.jl")

function fraction_within(rs, N, coords, hull)
    res = []
    vhull = VPolygon(hull)
    for r in rs
        m = zeros(N)
        for j in 1 : N
            circ, x, y, = sample_circle(centers(coords), r)
            m[j] = mean([[x[i], y[i]] ∈ vhull for i in eachindex(x)])
        end
        push!(res, (radius = r, frac = mean(m)))
    end
    return DataFrame(res)
end

makehull(coords) = convex_hull(map(r -> [r.x, r.y], eachrow(coords)))
mm(x) = minimum(x), maximum(x)

function sample_circle(coords, r)
    circ = circle_shape(coords[1], coords[2], r)
    θ = rand(Uniform(0, 2π), 100)
    x = coords[1] .+ r .* cos.(θ)
    y = coords[2] .+ r .* sin.(θ)
    return circ, x, y
end

function centers(df, n = 1)
    s = StatsBase.sample(1 : nrow(df), n)
    return df[s, :x], df[s, :y]
end

function circle_shape(h, k, r)
	θ = LinRange(0, 2 * π, 500)
	return h .+ r * sin.(θ), k .+ r * cos.(θ)
end

function diagplots(res)
    main = "Capacity at distance d: cap(d) = d * a2π\na: average fraction of circle with random
center and radius d that falls within Germany"
    p1 = scatter(res.radius, res.pot,
                 xlab = "Distance in km",
                 ylab = "Capacity: Distance * a2π",
                 title = main, label = "", size = (900, 600))
    p2 = scatter(res.radius, res.frac,
                 xlab = "Distance", ylab = "fraction within", label = "",
                 title = "Fraction of circle falling within Germany")
    return plot(p1, p2)
end

function barplots(df, bw = 20)
    df.distbin = floor.(df.dist ./ bw) .* bw  # Group into bins like [0,1), [1,2), etc.
    df2 = combine(groupby(df, :distbin), [:flows, :fpp] .=> sum)
    b1 = bar(df2.distbin, log.(df2.fpp_sum),
             xlab = "distance binned",
             ylab = "log(flows / a2pir)", label = "")

    b2 = bar(df2.distbin, log.(df2.flows_sum),
             xlab = "distance binned",
             ylab = "log(flows)", label = "")
    age = unique(df.agegroup)[1]
    year = unique(df.year)[1]
    return df2, plot(b1, b2, plot_title = "$age in $year")
end

function plothull(coords, add_pts = true)
    xmin, xmax = mm(coords[!, :x])
    ymin, ymax = mm(coords[!, :y])
    w = xmax - xmin
    h = ymax - ymin
    p = plot(VPolygon(hull), aspect_ratio = h / w)
    if add_pts; scatter!(coords.x, coords.y, markersize = 0.5, label = ""); end
    return p
end

function maincircle(df, districts)
    coords = select(unique(districts, :distcode),
                    [:xcoord => :x, :ycoord => :y])
    hull = makehull(coords)

    res = fraction_within(1:821, 10^3, coords, hull)
    res.pot = 2π .* res.radius .* res.frac
    diagplots(res)

    leftjoin!(df, res, on = [:dist => :radius])
    df.pot = Float64.(df.pot)
    df.fpp = df.flows ./ df.pot
    return res
end

# df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
# ## df = age(year(df, 2017), "30-50")
# districts = CSV.read("../data/districts.csv", DataFrame)

## res = maincircle(df, districts)

# groups = groupby(df, [:agegroup, :year])
# dfs, plots = barplots.(collect(groups), 25)

# phull = plothull(coords, true)
# for _ in 1:20
#     c, x, y = sample_circle(sample_coords(coords), rand(1:821))
#     plot!(phull, c, label = "")
# end
# phull
