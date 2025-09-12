using CSV, DataFrames, GeoStats, CairoMakie
include("./utils.jl")

outp = "/home/konstantin/paper/sections/texdocs/images/"

function eastwest(ags)
    east = r"11|12|13|15|16"
    if startswith(ags, east); return "East-Germany"; end
    if startswith(ags, "14"); return "Saxony"; end
    return "West-Germany"
end

uniqueN(x) = length(unique(x))

df = CSV.read("../../data/clean/correct.csv", DataFrame)
df = df[df.year .<= 2021, :]
east = r"11|12|13|14|15|16"
df.group = eastwest.(string.(df.ags_old));

dfp = combine(groupby(df, [:year, :group]), :ags_old => uniqueN => :N)
sax = year(df[df.group .== "Saxony", :], 2018).ags_new
dfp2 = combine(groupby(filter(:ags_new => n -> n ∈ sax, df), [:year, :ags_new, :name_new]),
              :ags_old => uniqueN => :N)
dfp2 = sort(year(dfp2, 1990), :N)
dfp2[dfp2.ags_new .== 14628, :name_new] .= "SSO*"

f = Figure(size = (400, 300), fontsize = 10);
ax = Axis(f[1, 1],
          xlabel = "Year",
          ylabel = "Number of districts",
          xgridvisible = false, ygridvisible = false)
addlines!(dfp, ax, dfp.group, :year, :N)
axislegend(ax, patchsize = (2, 10), framevisible = false, labelsize = 10,
           position = (.2, .8))

ax2 = Axis(f[1, 2],
           xlabel = "Number of districts in 1990",
           yticks = (1:nrow(dfp2), dfp2.name_new),
           xgridvisible = false, ygridvisible = false,
           )
scatter!(ax2, dfp2.N, 1:nrow(dfp2))
resize_to_layout!(f)
save(joinpath(outp, "adminchanges.pdf"), f)

dfp = df[df.conv_p .< 1, :]
f = Figure(size = (500, 300), fontsize = 10);
ax1 = Axis(f[1, 1], xlabel = "Share of population transferred (1990 - 1999)",
           ylabel = "Counts", limits = (0, 1, 0, 200))
plotgroups!(year(dfp, 1990.0, 1999.0), ax1, hist!, :group, :conv_p, nothing)
ax2 = Axis(f[1, 2], xlabel = "Share of population transferred (2000 - 2018)",
           ylabel = "Counts", limits = (0, 1, 0, 200))
plotgroups!(year(dfp, 2000.0, 2018.0), ax2, hist!, :group, :conv_p, nothing)
Legend(f[0, :], ax1, orientation = :horizontal)
save(joinpath(outp, "share_transferred.pdf"), f)

# f <- c(0, 100, 20, 20, 0, 100, 50, 10, 0)
# fm <- matrix(c(0, 20, 50, 100, 0, 10, 20, 100, 0), ncol = 3)
# o <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
# d <- c("A", "B", "C", "A", "B", "C", "A", "B", "C")
# flows <- data.table(origin = o, destination = d, flow = f, year = 2017)

# c <- c(1, 0, 0, 1, 0, 0, .2, 0, .8)
# cm <- matrix(c(1, 1, .2, 0, 0, 0, 0, 0, .8), ncol = 3)

# t(cm)
# fm
# sum(fm)
# (t(cm) %*% fm)
# t(cm) %*% fm %*% cm


# fm2 <- c(0, 10, 10, 100, 0,
#          100, 0, 20, 10, 0,
#          20, 100, 0, 10, 0,
#          10, 20, 100, 0, 0,
#          0, 0, 0, 0, 0)
# fm2 <- matrix(fm2, ncol = 5, byrow = TRUE)
# cm2 <- c(1, 0, 0, 0, 0, .2, .8, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
# cm2 <- matrix(cm2, ncol = 5, byrow = TRUE)

# new <- t(cm2) %*% fm2 %*% cm2
