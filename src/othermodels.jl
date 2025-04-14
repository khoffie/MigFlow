function distonly(data::NamedTuple)
    flows = data.flows
    fromdist = data.fromdist
    todist = data.todist
    fp = data.frompop
    tp = data.topop
    di = data.dist
    ds = data.distscale
    m = maximum(di)
    @model function model(flows, fp, tp, di, ds, m)
        a ~ Normal(-8, 1)
        c ~ Gamma(15, .2) ## we want mean ~ 3 but with less variance than Gamma(3, 1)
        l ~ Gamma(10, 1)
        d0 ~ Gamma(10, 1)

        c = c / 10 # to offest smaller variance of c
        l = l / 100
        d0 = d0 / 100
        di = di ./ ds
        tp = tp ./ 153000 # median pop

        pop = log.(tp)
        dist = log.(l .* 1 .+ (1 - l) ./ ((di .+ d0).^c))

        preds = fp .* logistic.(a .+ pop .+ dist)
        flows ~ arraydist(Poisson.(preds))
        return preds
    end

    mdl = model(flows, fp, tp, di, ds, m)
    lb = [-20, 10, 0, 1]
    ub = [0, 100, 99, 100]
    return (; mdl, lb, ub)
end

function t(x, c, m = 100.0, e = 0.1)
    B = (c + 1) / ((m + e)^(c+1) - e^(c+1))
    return B * (x + e)^c
end

function gravity(data::NamedTuple)
    flows = data.flows
    fp = data.frompop
    tp = data.topop
    di = data.dist
    ds = data.distscale

    @model function model(flows, fp, tp, di, ds = 100)
        a ~ Normal(-8, 1)
        c ~ Gamma(15, .2) ## we want mean ~ 3 but with less variance than Gamma(3, 1)
        s ~ Gamma(1, 1)
        t ~ Gamma(1, 1)

        c = c / 10 # to offest smaller variance of c
        di = di ./ ds
        tp = tp ./ 153000 # median pop

        pop = log.(tp.^t)
        dist = log.(1 ./ (di).^c)
        preds = fp.^s .* logistic.(a .+ pop .+ dist)

        flows ~ arraydist(Poisson.(preds))
        return preds
    end

    mdl = model(flows, fp, tp, di, ds)
    lb = [-20, 10, 0, 0]
    ub = [0, 100, 5, 5]
    return (; mdl, lb, ub)
end

function gravity_nopopc(data::NamedTuple)
    flows = data.flows
    fp = data.frompop
    tp = data.topop
    di = data.dist
    ds = data.distscale

    @model function model(flows, fp, tp, di, ds = 100)
        a ~ Normal(-8, 1)
        c ~ Gamma(15, .2) ## we want mean ~ 3 but with less variance than Gamma(3, 1)

        c = c / 10 # to offest smaller variance of c
        di = di ./ ds
        tp = tp ./ 153000 # median pop

        pop = log.(tp)
        dist = log.(1 ./ (di).^c)
        preds = fp .* logistic.(a .+ pop .+ dist)

        flows ~ arraydist(Poisson.(preds))
        return preds
    end

    mdl = model(flows, fp, tp, di, ds)
    lb = [-20, 10]
    ub = [0, 100]
    return (; mdl, lb, ub)
end

# @model function linear(flows, fp, tp, di, ds = 100)
#     a ~ Normal(-8, 1)
#     b ~ Gamma(1, 5)
#     c ~ Gamma(15, .2) ## we want mean ~ 3 but with less variance than Gamma(3, 1)

#     c = c / 10 # to offest smaller variance of c
#     di = di ./ ds
#     tp = tp ./ 153000 # median pop

#     pop = log.(tp)
#     dist = log.(1 ./ (b .+ c .* di))
#     preds = fp .* logistic.(a .+ pop .+ dist)

#     flows ~ arraydist(Poisson.(preds))
#     return preds
# end
