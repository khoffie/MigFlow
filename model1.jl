@model function migration1old(flows,fromdist,todist,frompop,topop,distance,agegroup,Ndist,meddist)
    a ~ Gamma(5.0, 1.0/4.0)
    b ~ Gamma(3.0, 10.0/2.0)
    c ~ Gamma(5.0, 1.0/4.0)
    d0 ~ Gamma(5.0, 0.10/4.0)
    preds = frompop .* topop .* a ./ 1000.0 .* (1.0 .+ b ./ (distance ./ meddist .+ d0).^c)
    flows ~ arraydist([Poisson(p) for p in preds])
    return preds
end

@model function migration1new(flows, fromdist, todist, frompop, topop, distance,
                              agegroup, Nages, Ndist, meddist)
    a ~ filldist(Gamma(5.0, 10.0/4.0), Nages)
    b ~ filldist(Gamma(3.0, 1.0/2.0), Nages)
    c ~ filldist(Gamma(5.0, 1.0/4.0), Nages)
    d0 ~ filldist(Gamma(5.0, 0.10/4.0), Nages)

    preds = [frompop[i] * topop[i] * a[agegroup[i]] / 1000.0 * (1.0 + b[agegroup[i]] / (distance[i] / meddist + d0[agegroup[i]])^c[agegroup[i]]) for i in 1:length(flows)]

    if typeof(a[1]) != Float64
        @show a[1].value,b[1].value,c[1].value,d0[1].value
    end
    if any(isnan,preds)
        println("NaN in predictions")
    end
    flows ~ arraydist([Poisson(p) for p in preds])
    return preds
end

# sample the above model for 600 samples on 3 parallel chains with NUTS after doing 500 warmups and targeting 80% acceptance
# inis = [3.54; 0.25; 3.85; 0.13]
# fit1 = sample(migration1(dt2.flows,
#                          levelcode.(dt2.fromdist),
#                          levelcode.(dt2.todist),
#                          dt2.frompop_ger,
#                          dt2.topop,
#                          dt2.distance,
# ##                         ourdat.gdpcfrom, ourdat.gdpcto,
#                          levelcode.(dt2.agegroup),
#                          Ndist,
#                          meddist),
#     NUTS(500,.8), MCMCThreads(), 600, 3; init_params = Iterators.repeated(inis))

# plot(fit1)
