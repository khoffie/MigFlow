using CSV,DataFrames, DataFramesMeta, Distributions, Random,LinearAlgebra, CategoricalArrays, StatsPlots


function simdataset!(rng, pairsdata,simrate)

    pairsdata.flow = [rand(rng,Poisson(simrate(r))) for r in eachrow(pairsdata)]
        
end

function abcrate(row)
    a = 11.36
    b = 3.56
    c = 1.92
    d0 = 0.03

    row.frompop * row.topop * a/1000 * (1 + b/(row.distance/row.meddist + d0)^c) * row.todesir / row.fromdesir
end

function makeflowstemplate(popfile,centroidfile,outfile)

    munipop = CSV.read(popfile,DataFrame)
    munipopsum = @by(munipop, :district,:totpop = sum(:pop) ./ 1000)
    munipopsumage = @by(munipop, [:district,:age_group],:totpop = sum(:pop) ./ 1000)
    println("There are $(nrow(munipopsum)) districts and $(nrow(munipopsumage)) district/age pairs")

    municent = CSV.read(centroidfile,DataFrame)
    distcent = @chain leftjoin(munipopsum,municent,on = :district) begin
        @by(:district, :xcent = sum(:totpop .* :centroid_x), :ycent = sum(:totpop .* :centroid_y))
    end
    distcent = leftjoin(distcent,munipopsum,on=:district)
    distcent.xcent = distcent.xcent ./ distcent.totpop
    distcent.ycent = distcent.ycent ./ distcent.totpop
    select!(distcent,Not(:totpop))
    templflow = leftjoin(munipopsumage,distcent; on = :district, makeunique=true)

    templflow = @subset(crossjoin(templflow,templflow; makeunique=true), :district .!= :district_1,:age_group .== :age_group_1)
    @show templflow[1:10,:]
    rename!(templflow,Dict(:district => :fromdist, :district_1 => :todist,:totpop => :frompop, :totpop_1 => :topop))
    @transform!(templflow,@byrow :distance = norm([:xcent,:ycent] - [:xcent_1,:ycent_1])/1e6)
    @select!(templflow,:fromdist,:todist,:distance,:frompop,:topop,:age_group)
    CSV.write(outfile,templflow)

end

if false
    makeflowstemplate("data/munis_pop.csv","data/munis_centroid.csv","data/flowtemplate.csv")

end

function simdatafromtemplate(rng,templfile)

    flowtemp = CSV.read(templfile,DataFrame)
    flowtemp.fromdist = categorical(flowtemp.fromdist)
    flowtemp.todist = categorical(flowtemp.todist)
    desir = rand(rng,Gamma(5.0,1.0/4.0),length(levels(flowtemp.todist)))
    flowtemp.fromdesir = desir[levelcode.(flowtemp.fromdist)]
    flowtemp.todesir = desir[levelcode.(flowtemp.todist)]
    meddist = median_distance()
    flowtemp.meddist = fill(meddist,nrow(flowtemp))
    simdataset!(rng,flowtemp,abcrate)
    (flows = flowtemp, desirability = desir)
end

if false

    simdata = simdatafromtemplate(Xoshiro(20240725),"data/flowtemplate.csv")


end