using CategoricalArrays

ages = ["below18", "18-25", "25-30","30-50", "50-65", "above65"]
flows = read_flows("../data/FlowDataGermans.csv", 1.0, ages, 2000:2017)
flows = flows[(flows.agegroup .== "30-50") .& (flows.year .== 2017), :]

districts = CSV.read("../data/districts.csv", DataFrame)

add_lrd(districts)
flows = joinlrd(flows, districts)


mdat = gen_mdat(flows, districts; distscale = 100.0, ndc = 1, ngc = 1)
out = estimate(distonly, mdat)

df = DataFrame(flows = mdat.flows, preds = out.preds,
               fromdens = mdat.fromdens, todens = mdat.todens)
df.res = df.flows ./ df.preds
mind, maxd = minmax(df.fromdens)
ng = 10
df.odens = cut(df.fromdens, ng; labels = ["Q$i" for i in 1:ng])
df.tdens = cut(df.todens, ng; labels = ["Q$i" for i in 1:ng])
hm = combine(groupby(df, [:odens, :tdens]), :res => mean)
hmm = unstack(hm, :odens, :tdens, :res_mean)
heatmap(Matrix(hmm[:, 2:end]))


Matrix{Float64}(hmm)



heatmap(hmm)


resm = unstack(df, :fromdens, :todens, :res, combine = last)

heatmap(Matrix(resm))

?unstack


density(df.res)
df
s = length(unique(df.fromdens))

reshape(df.res, (s, s))



function plotgeocheby(ncoefs)
    geo, p = evalgeocheby(rand(ncoefs), unique(districts, :distcode), true)
    p1 = scatter(geo.xcoord, geo.geo, xlab = "xcoord")
    p2 = scatter(geo.ycoord, geo.geo, xlab = "ycoord")
    p3 = scatter(geo.ycoord .* geo.xcoord, geo.geo, xlab = "xcoord * ycoord")
    plot(p, p1, p2, p3, plot_title = "Geocheby with $ncoefs coefs")
end



mu = fill(0.0, length(vals.funval))
sigma = fill(3.0, length(vals.funval))
## logpdf returns the log of density of MvN at all 401 districts
Turing.@addlogprob!(logpdf(MvNormal(mu, sigma), desvals))
## tamp down the corners of the map
Turing.@addlogprob!(logpdf(MvNormal(fill(0.0, 4), fill(0.125, 4)),
                           [desfun(xmin, ymin),
                            desfun(xmin, ymax),
                            desfun(xmax, ymin),
                            desfun(xmax, ymax)]))
