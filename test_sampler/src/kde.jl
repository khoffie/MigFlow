using Distributions, KernelDensity, StatsPlots

dist = SkewNormal(0, 10, .105)
dens = kde([1, 10, 100], dist)
plot(dens)
dist = SkewNormal(0, 10, .106)
kde([1, 10, 100], dist)
