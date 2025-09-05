
######################################################################
##################### Data Issues ####################################
df1 = destination(origin(df, [3159]), [5978])
combine(groupby(df1, :year), :flows => sum)

sort(outflux(origin(df, 3159), [:year]), :outflux)
sort(outflux(destination(origin(df, 3159), 5978), [:year]), :outflux)

######################################################################
############################# GLM ####################################

df17 = combine(groupby(year(df, 2017), [:fromdist, :todist]),
               [:flows => sum => :flows,
                :topop => first => :topop,
                :dist => first => :dist,
               :frompop => sum => :frompop])
leftjoin!(df17, year(di, 2017)[!, [:distcode, :area]],
          on = :todist => :distcode)
df17.toarea = Float64.(df17.area)
select!(df17, Not(:area))
leftjoin!(df17, year(di, 2017)[!, [:distcode, :area]],
          on = :fromdist => :distcode)
df17.fromarea = Float64.(df17.area)
select!(df17, Not(:area))

sum(unique(df17, :topop).topop)
sum(unique(df17, :frompop).frompop)

m = glm(@formula(flows ~ log(dist)), df17, Poisson(), LogLink())
m = glm(@formula(flows ~ log(frompop) + log(topop) + log(dist)), df17, Poisson(), LogLink())
1 - var(df17.flows .- predict(m)) / var(df17.flows)
