using CSV, DataFrames, CairoMakie, StatsBase
outp = "/home/konstantin/paper/sections/texdocs/images/"

distcode(m) = length(m) == 7 ? SubString(m, 1, 4) : SubString(m, 1, 5)
statecode(m) = length(m) == 4 ? SubString(m, 1, 1) : SubString(m, 1, 2)

district(df, d) = filter(row -> row.distcode ∈ d, df)
calcdensity(md, ma) = mean(md, Weights(ma))
gb = DataFrames.groupby

function readinkar()
    f = "/home/konstantin/Backup/workstation/home/Diss/inst/extdata/clean/inkar/inkar_2021.csv"
    ink = CSV.read(f, DataFrame)
    ink = ink[ink.Raumbezug .== "Gemeinden", :]
    return ink[ink.Zeitbezug .== "2017", :]
end

function munidf(ink)
    df = unstack(ink[ink.ID .== 425 .|| ink.ID .== 426, :], :Kennziffer, :ID, :Wert)
    rename!(df, "425" => :pop)
    rename!(df, "426" => :area)
    rename!(df, :Kennziffer => :muni)
    df.distcode = distcode.(df.muni)
    df.pop  = parse.(Float64, replace.(df.pop, ',' => '.'))
    df.area = parse.(Float64, replace.(df.area, ',' => '.'))
    df.muni = parse.(Int, df.muni)
    df.distcode = parse.(Int, df.distcode)
    df.density = df.pop ./ df.area * 100
    return filter(row -> row.pop > 0, df)
end

function districtdf(df)
    means = combine(gb(df, :distcode), [:density, :area] => calcdensity => :mean)
    stds = combine(gb(df, :distcode), :density => std => :sd)
    stds.sd = map(x -> isnan(x) ? 0.0 : x, stds.sd)
    df2 = leftjoin(means, stds, on = :distcode)
    df2.rel = df2.mean ./ df2.sd
    df2.pop = combine(gb(df, :distcode), :pop => sum => :pop).pop
    return df2
end

function violinplot(munis, districts, distcodes)
    m = "Distribution of population density of\nmunicipalities within districts"
    mu = district(munis, distcodes)
    di = district(districts, distcodes)
    violin(string.(mu.distcode), mu.density, label = "",
           ylab = "Population density", xlab = "",
           title = m, xticks = false)
    dotplot!(string.(mu.distcode), mu.density,
             marker = (:black, stroke(0), log.(mu.pop) ./ 3), label = "")
    scatter!(string.(di.distcode), di.mean,
            marker = (:red, stroke(0), log.(di.pop) ./ 3),
            label = "")
    savefig("../../mig-paper/images/density_munis.svg")
end

df = munidf(readinkar())
dis = districtdf(df)
samples = StatsBase.sample(dis[dis.sd .> 0.0, :].distcode, 10)
mu = district(df, samples)
di = district(dis, samples)
mu.group = groupindices(groupby(mu, :distcode))
di.group = groupindices(groupby(di, :distcode))

f = Figure(size = (400, 200), fontsize = 10);
ax = Axis(f[1, 1], xlabel = "", ylabel = "Population Density",
          title = "Distribution of population density of\nmunicipalities within districts",
          xgridvisible = false, ygridvisible = false,
          xticksvisible = false )
ylims!(ax, 0, nothing)
violin!(ax, mu.group, mu.density, width = 2.0)
xjitter = mu.group .+ 0.2 .* (rand(length(mu.group)) .- 0.5)
scatter!(ax, xjitter, mu.density, color = :black, markersize = log.(mu.pop) ./ 3)
scatter!(ax, Int.(di.group), di.mean, color = :red, markersize = log.(di.pop) ./ 3)
hidexdecorations!(ax)
save(joinpath(outp, "density_munis.pdf"), f)
