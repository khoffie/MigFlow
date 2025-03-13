function distonly(data::NamedTuple)
    flows = data.flows
    fp = data.frompop
    tp = data.topop
    di = data.dist
    ds = data.distscale

    @model function model(flows, fp, tp, di, ds)
        a ~ Normal(-8, 1)
        c ~ Gamma(15, .2) ## we want mean ~ 3 but with less variance than Gamma(3, 1)
        l ~ Gamma(10, 1)
        d0 ~ Gamma(10, 1)
        t ~ Gamma(1, 1)
        e ~ Normal(1, 1)

        c = c / 10 # to offest smaller variance of c
        l = l / 100
        d0 = d0 / 100
        di = di ./ ds
        tp = tp ./ 153000 # median pop

        pop = log.((tp .+ e).^t)
        dist = log.(l .+ (1 - l) ./ (di .+ d0).^c)
        preds = fp .* logistic.(a .+ pop .+ dist)

        flows ~ arraydist(Poisson.(preds))
        return preds
    end

    mdl = model(flows, fp, tp, di, ds)
    lb = [-20, 10, 0, 1, 0, 0]
    ub = [0, 100, 99, 100, 5, 10]
    return (; mdl, lb, ub)
end

function distonlynoe(data::NamedTuple)
    flows = data.flows
    fp = data.frompop
    tp = data.topop
    di = data.dist
    ds = data.distscale

    @model function model(flows, fp, tp, di, ds)
        a ~ Normal(-8, 1)
        c ~ Gamma(15, .2) ## we want mean ~ 3 but with less variance than Gamma(3, 1)
        l ~ Gamma(10, 1)
        d0 ~ Gamma(10, 1)
        t ~ Gamma(1, 1)

        c = c / 10 # to offest smaller variance of c
        l = l / 100
        d0 = d0 / 100
        di = di ./ ds
        tp = tp ./ 153000 # median pop

        pop = log.(tp.^t)
        dist = log.(l .+ (1 - l) ./ (di .+ d0).^c)
        preds = fp .* logistic.(a .+ pop .+ dist)

        flows ~ arraydist(Poisson.(preds))
        return preds
    end

    mdl = model(flows, fp, tp, di, ds)
    lb = [-20, 10, 0, 1, 0]
    ub = [0, 100, 99, 100, 5]
    return (; mdl, lb, ub)
end

function distdens(data::NamedTuple)
    flows = data.flows
    fp = data.frompop
    tp = data.topop
    fd = data.fromdens
    td = data.todens
    di = data.dist
    ds = data.distscale
    ndc = data.ndenscoefs
    dmin = minimum(data.fromdens)
    dmax = maximum(data.fromdens)

    @model function model(flows, fp, tp, di, fd, td,
                          ndenscoefs, densmin, densmax, ds)
        a ~ Normal(-8, 1)
        c ~ Gamma(15, .2) ## we want mean ~ 3 but with less variance than Gamma(3, 1)
        l ~ Gamma(10, 1)
        d0 ~ Gamma(10, 1)
        kd ~ MvNormal(zeros(ndenscoefs), fill(10.0, ndc))
        c = c / 10 # to offest smaller variance of c
        l = l / 100
        d0 = d0 / 100
        kd = kd / 10
        kd[1] = 0.0 # cheby intercept, ensure heatmap is not elevated
        di = di ./ ds
        tp = tp ./ 153000 # median pop

        pop = log.(tp)
        dist = log.(l .+ (1 - l) ./ (di .+ d0).^c)
        densitycheby = defdensitycheby(kd, dmin, dmax)
        density = densitycheby.(fd, td)

        preds = fp .* logistic.(a .+ pop .+ dist .+ density)

        flows ~ arraydist(Poisson.(preds))
        return preds
    end

    mdl = model(flows, fp, tp, di, fd, td, ndc, dmin, dmax, ds)
    lb = [-20, 10, 0, 1, fill(-10, ndc)...]
    ub = [0, 100, 99, 100, fill(10, ndc)...]
    return (; mdl, lb, ub)
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
