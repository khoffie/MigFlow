using CSV, DataFrames, Revise, CairoMakie, StatsBase, GeoStats

outp = "/home/konstantin/paper/sections/texdocs/images/"
includet("utils.jl")

function addlines!(df, ax, group, xcol, ycol)
    for g in unique(group)
        foo = df[group .== g, :]
        lines!(ax, foo[!, xcol], foo[!, ycol], label = a)
    end
end


df = CSV.read("../../data/FlowDataGermans.csv", DataFrame);
df2 = outflux(df, sum, [:todist, :year])
out = outflux(df2, sum, [:year])
districts = CSV.read("../../data/districts.csv", DataFrame);
correct = CSV.read("../../data/clean/correct.csv", DataFrame)
ages = unique(df.agegroup)
allyears = unique(df.year)

## In the early 2000s are unreasonably high flows. What is going on?
## To get a quick idea we will count the flows that are 5x higher than
## the average between 2006-2008

function toolarge(df, flow_th, avg_th)
    dub = df2[df2.year .<= 2005, :]
    ave = outflux(year(df2, 2006.0, 2008.0), mean, [:todist], "average")
    leftjoin!(dub, ave, on = [:fromdist, :todist])
    dub.rel = dub.flows ./ dub.average
    dub.toolarge = dub.rel .>= flow_th
    dub = dub[dub.average .>= avg_th, :]
    return dub, sort(combine(groupby(dub, [:fromdist, :year]), :toolarge => sum => :N), :N)
end

dub, dflarge = toolarge(df2, 5.0, 5.0)

## https://www.bva.bund.de/DE/Services/Buerger/Migration-Integration/Spaetaussiedler/04_Informationen/Grenzdurchgangslager/Grenzdurchgangslager_node.html?utm_source=chatgpt.com
## Here 9 centers are mentioned, which existed until october 2000:
## Friedland, Nürnberg, Osnabrück, Bramsche, Hamm, Empfingen, Dranse
## und Schönberg-Holm. (the source says 9 but only lists 8)

##Friedland in Göttingen, Bramsche in district Osnabrück, Empfingen in
## Freudenstadt, Dranse in Ostprignitz-Ruppin, Schönberg-Holm in Plön

late_names = ["Göttingen", "Nürnberg" , "Osnabrück Kreis", "Hamm", "Freudenstadt", "Plön"]
late = DataFrame(name = late_names)
leftjoin!(late, year(districts, 2017)[!, [:distcode, :name]], on = :name)
late.distcode = Int.(late.distcode)

f = Figure(size = (600, 400), fontsize = 10);
ags = reshape(late.distcode, 2, 3)
for i in 1:size(ags, 1)
    for j in 1:size(ags, 2)
        c = ags[i, j]
        ax = Axis(f[i, j],
                  xlabel = "Year",
                  ylabel = "Outflux (in thousands)",
                  title = code(late, c).name[1],
                  xgridvisible = false, ygridvisible = false)
        df = origin(out, c)
        lines!(ax, df.year, df.flows ./ 1000)
    end
end
save(joinpath(outp, "late.pdf"), f)
