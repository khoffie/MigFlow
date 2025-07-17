using Plots, StatsBase, ImageFiltering, Random, Distributions

function genmat(N; cx = 500, cy = 500, pmove = .1)
    D = 1000
    m = zeros(Float64, (D, D))
    for i in 1:N
        S = rand(10 : 50)
        xidx = rand(1 : D - S)
        yidx = rand(1 : D - S)
        dist = sqrt((xidx - cx)^2 + (yidx - cy)^2)
        μ = (1 + sqrt(dist))^(-1)
        m[xidx : xidx + S, yidx : yidx + S] .= μ ## .+ rand(Gamma(1, 1))
    end
    return m ./ sum(m) .* pmove
end

heat(m, max) = heatmap(m, clim = (0.0, max), colorbar = false, xaxis = false, yaxis = false)

function plotheatmaps()
    m1 = genmat(100; pmove = .1)
    m2 = genmat(100; cx = 800, cy = 700, pmove = .12)
    maxval = max(maximum(m1), maximum(m2))

    m1blurred = imfilter(m1, Kernel.gaussian(200))
    m2blurred = imfilter(m2, Kernel.gaussian(200))
    maxvalblurred = max(maximum(m1blurred), maximum(m2blurred))

    U = ones(1000, 1000)
    plot(heat(m1, maxval), heat(m1blurred, maxvalblurred), heat(sum(m1) .* U, max(sum(m1), sum(m2))),
         heat(m2, maxval), heat(m2blurred, maxvalblurred), heat(sum(m2) .* U, max(sum(m1), sum(m2)))
         )
end
