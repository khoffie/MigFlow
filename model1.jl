@model function migration1(flows,fromdist,todist,frompop,topop,distance,agegroup,Ndist,meddist)
    a ~ Gamma(5.0, 1.0/4.0)
    b ~ Gamma(3.0, 10.0/2.0)
    c ~ Gamma(5.0, 1.0/4.0)
    d0 ~ Gamma(5.0, 0.10/4.0)
    preds = frompop .* topop .* a ./ 1000.0 .* (1.0 .+ b ./ (distance ./ meddist .+ d0).^c)
    flows ~ arraydist([Poisson(p) for p in preds])
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
