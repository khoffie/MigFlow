using Pkg
Pkg.activate(".")

using CSV,DataFrames,DataFramesMeta,Random,Distributions,StatsBase, StatsFuns,
    StaticArrays,LinearAlgebra,StatsPlots,GLM, CategoricalArrays


munipop = CSV.read("./data/munis_pop.csv",DataFrame)

munidist = @by(munipop,:municipality,:district = :district[1])

munitotpop = @by(munipop, [:municipality,:age_group], :totpop = sum(:pop))
muniloc = CSV.read("./data/munis_centroid.csv",DataFrame)
muniloc = leftjoin(muniloc,munidist; on=:municipality)

function simmoves(rn,n,munipop,locs,distmap)

    dat = DataFrame()

    for age in unique(munipop.age_group)
        muni = @subset(munipop,:age_group .== age)
        pops = Dict(muni.municipality .=> muni.pop)
        ws = FrequencyWeights(muni.pop)
        munis = muni.municipality
        from = sample(rn,munis,ws)
        frompop = pops[from]
        to = sample(rn,munis,ws)
        topop = pops[to]

        for i in 1:n
            newfrom = sample(rn,munis,ws)
            newfrompop = pops[newfrom]
            newto = sample(rn,munis,ws)
            newtopop = pops[newto]
            if newtopop > newfrompop || rand(rn) < newtopop / newfrompop
                to,from = newto,newfrom
                topop,frompop = newtopop,newfrompop
            end
            iscaptured = distmap[from] != distmap[to]
            push!(dat,(from = from,to = to, dist = norm(locs[from]-locs[to]),captured=iscaptured,agegroup = age, fromdist = distmap[from], todist=distmap[to]))
        end
    end
    dat
end

locs = Dict(muniloc.municipality .=> [SVector(r.centroid_x,r.centroid_y) for r in eachrow(muniloc)])
distmap = Dict(muniloc.municipality .=> muniloc.district)

r = Xoshiro(20240709)

movdat = simmoves(r,100000,munipop,locs,distmap)

p = @df movdat scatter(:dist ./ 1000.0,:captured,xlab="Distance Moved (km)",ylab="Captured Move? (y/n)",
    title="Probability To Detect by Distance",xlim=(0,100))


movdat.agegroupcat = categorical(movdat.agegroup)
movdat.fromdistcat = categorical(movdat.fromdist)
movdat.todistcat = categorical(movdat.todist)

mod1 = glm(@formula(captured ~  agegroupcat + dist/100000.0 & agegroupcat),movdat,Bernoulli())
mod2 = glm(@formula(captured ~ fromdistcat + agegroupcat + dist/100000.0 & agegroupcat),movdat,Bernoulli())
@show(mod1)
@show(mod2)


for agegrp in levels(movdat.agegroup)

    predfor = DataFrame(dist=collect(0.0:100:100000),agegroupcat=categorical(fill(agegrp,1001)))
    predfor.preds = predict(mod1,predfor)
    p = @df predfor plot!(:dist ./ 1000.0,:preds;label=agegrp)
end
display(p)



## By district, show the distribution of detected move size:
p=density()
for i in unique(movdat.fromdist)
    det = @subset(movdat[movdat.fromdist .== i .&& movdat.captured,:])
    global p = density!(det.dist./1000.0; title="Detected Moves (each district by color)",
        alpha=.4,xlab="Distance (km)",xlim=(0,1000),legend=false)
end
display(p)



density(movdat[.! movdat.captured,:dist]./1000.0; title="Undetected moves (all districts)",xlab="distance (km)")