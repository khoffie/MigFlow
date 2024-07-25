function desir(d1,d2,d3,d4,d5,d6,age,from,to)
    if age == 1
        d1[to]/d1[from]
    elseif age == 2
        d2[to]/d2[from]
    elseif age == 3
        d3[to]/d3[from]
    elseif age == 4
        d4[to]/d4[from]
    elseif age == 5
        d5[to]/d5[from]
    elseif age == 6
        d6[to]/d6[from]                
    end
end

@model function migration2(flows,fromdist,todist,frompop,topop,distance,agegroup,Nages,Ndist,meddist, netactual)
    a ~ Gamma(5.0,10.0/4.0)
    b ~ Gamma(3.0,1.0/2.0)
    c ~ Gamma(5.0,1.0/4.0)
    d0 ~ Gamma(5.0,0.10/4.0)
    neterr ~ Gamma(3.0,0.1/2.0)

    desir1 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir2 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir3 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir4 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir5 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir6 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desires = [desir(desir1,desir2,desir3,desir4,desir5,desir6,
                        agegroup[i],fromdist[i],todist[i]) for i in 1:length(flows)]
    preds = frompop .* topop .* a ./ 1000.0 .* (1.0 .+ b ./ (distance ./ meddist .+ d0).^c) .* desires

    if typeof(a) != Float64
        @show a.value,b.value,c.value,d0.value,neterr.value
    end
    if any(isnan,preds) println("NaN in predictions") end
    netflows = calcnet(preds,fromdist,todist,agegroup,Nages,Ndist)

    flows ~ arraydist([Poisson(p) for p in preds])
    netactual ~ arraydist(Normal.(netflows,neterr .* abs.(netflows)))
    return((preds,netflows))
end

function calcnet(flows,fromdist,todist,agegrp,Nages,Ndist)
    netflows = zeros(typeof(flows[1]),(Ndist,Nages))
    for n in 1:length(flows)
        netflows[fromdist[n],agegrp[n]] -= flows[n]
        netflows[todist[n],agegrp[n]] += flows[n]
    end
    netflows
end
