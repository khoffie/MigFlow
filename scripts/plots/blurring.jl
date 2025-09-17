using ImageFiltering, Plots, LaTeXStrings
outp = "/home/konstantin/paper/sections/texdocs/images/"

d = 1:10
K(x, c, b) = exp(- (x - c)^2 / 2b^2)
w1 = K.(d, 5, 1)
w2 = K.(d, 5, 2)

f = Figure(size = (400, 300), fontsize = 10);
m = L"\text{Weights for Gaussian Kernel} \quad K(x, b) = \exp(- \frac{(x - 5)^2}{2b^2})"
ax = Axis(f[1, 1], title = m, ylabel = "Relative Contribution")
Makie.barplot!(ax, 1:10, w1 ./ sum(w1[w1 .> -Inf]), alpha = .1,
               label = "b = 1")
Makie.barplot!(ax, 1:10, w2 ./ sum(w2[w2 .> -Inf]), alpha = .1,
               label = "b = 10")
axislegend(ax)
save(joinpath(outp, "weights.pdf"), f)

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

heat(m, max) = Plots.heatmap(m, clim = (0.0, max), colorbar = false, xaxis = false, yaxis = false)

function plotheatmaps()
    m1 = genmat(100; pmove = .1)
    m2 = genmat(100; cx = 800, cy = 700, pmove = .12)
    maxval = max(maximum(m1), maximum(m2))

    m1blurred = imfilter(m1, Kernel.gaussian(200))
    m2blurred = imfilter(m2, Kernel.gaussian(200))
    maxvalblurred = max(maximum(m1blurred), maximum(m2blurred))

    U = ones(1000, 1000)
    p = Plots.plot(heat(m1, maxval), heat(m1blurred, maxvalblurred), heat(sum(m1) .* U, max(sum(m1), sum(m2))),
             heat(m2, maxval), heat(m2blurred, maxvalblurred), heat(sum(m2) .* U, max(sum(m1), sum(m2)))
             )
    return p
end

savefig(plotheatmaps(), joinpath(outp, "heat.pdf"))
