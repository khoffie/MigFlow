@model function migration1(flows, fromdist, todist, frompop, topop, distance,
                              agegroup, Nages, Ndist, meddist)
    a ~ filldist(Gamma(5.0, 10.0/4.0), Nages)
    b ~ filldist(Gamma(3.0, 1.0/2.0), Nages)
    c ~ filldist(Gamma(5.0, 1.0/4.0), Nages)
    d0 ~ filldist(Gamma(5.0, 0.10/4.0), Nages)

    preds = [frompop[i] * topop[i] * a[agegroup[i]] / 1000.0 * (1.0 + b[agegroup[i]] / (distance[i] / meddist + d0[agegroup[i]])^c[agegroup[i]]) for i in 1:length(flows)]

    if typeof(a[1]) != Float64
        @printf "a1 = %.2f, b1 = %.2f, c1 = %.2f, d0_1 = %.2f\n"  a[1].value b[1].value c[1].value d0[1].value
    end
    if any(isnan,preds)
        println("NaN in predictions")
    end
    flows ~ arraydist([Poisson(p) for p in preds])
    return preds
end
