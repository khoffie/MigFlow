function full2(data::NamedTuple; normalize = true)
    flows = data.flows
    fromdist = data.fromdist
    todist = data.todist
    fp = data.frompop
    tp = log.(data.topop ./ 153000) # median pop
    fd = data.fromdens
    td = data.todens
    di = data.dist ./ data.distscale
    ndc = data.ndc
    ngc = data.ngc
    dmin = data.dmin
    dmax = data.dmax
    distcode = data.distcode
    xcoord = data.xcoord
    ycoord = data.ycoord
    nfrom = length(unique(fromdist))
    N = length(flows)
    xmin, xmax = mm(xcoord)
    ymin, ymax = mm(ycoord)

    @model function model(flows, fromdist, todist, fp, tp, di, fd, td, N, nfrom,
                          ndc, ngc, dmin, dmax, xcoord, ycoord, distcode,
                          xmin, xmax, ymin, ymax)
        c ~ Gamma(15, .2); c = c / 10
        l ~ Gamma(10, 1); l = l / 100
        d0 ~ Gamma(10, 1); d0 = d0 / 100
        kd ~ MvNormal(zeros(ndc), fill(10.0, ndc)); kd = kd / 100
        kg ~ MvNormal(zeros(ngc), fill(10.0, ngc)); kg = kg / 100
        kd[1] = 0.0 # cheby intercept, ensure heatmap is not elevated
        kg[1] = 0.0

        T = eltype(c)  # to get dual data type for AD
        att   = Vector{T}(undef, N)
        denom = zeros(T, nfrom)
        preds = Vector{T}(undef, N)

        dist = log.(l .+ (1 - l) ./ (di .+ d0).^c)
        density = defdensitycheby(kd, dmin, dmax).(fd, td)
        geo = defgeocheby(kg, xmin, xmax, ymin, ymax).(xcoord, ycoord)

        @inbounds for i in 1:N
            att[i] = exp(tp[i] + dist[i] + density[i] +
                (geo[todist[i]] - geo[fromdist[i]]))
            denom[fromdist[i]] += att[i]
        end
        @inbounds for i in 1:N
                preds[i] = fp[i] * att[i] / denom[fromdist[i]]
        end
        flows .~ Poisson.(preds)
        return preds
    end

    mdl = model(flows, fromdist, todist, fp, tp, di, fd, td, N, nfrom,
                ndc, ngc, dmin, dmax, xcoord, ycoord, distcode,
                xmin, xmax, ymin, ymax)
    lb = [10, 0, 1, fill(-100, ndc)..., fill(-100, ngc)...]
    ub = [100, 99, 100, fill(100, ndc)..., fill(100, ngc)...]
    return (; mdl, lb, ub)
end
