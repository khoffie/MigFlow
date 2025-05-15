function norm(data::NamedTuple; normalize = true)
    flows    = data.flows
    fromdist = data.fromdist
    fp       = data.frompop
    tp       = data.topop ./ 153000 # median topop
    di       = data.dist  ./ data.distscale
    nfrom    = length(unique(fromdist))
    N        = length(flows)
    fpt      = data.fpt
    radius   = data.radius
    @model function model(flows, fromdist, fp, tp, di, nfrom, N,
                          fpt, radius, normalize)
        # Priors
        a      ~ Normal(-5, 1)
        b  ~ Gamma(1, 1);     ## b  = b_raw / 100
        c_raw  ~ Gamma(15, 0.2);  c  = c_raw / 10
        l_raw  ~ Gamma(10, 1.0);  l  = l_raw / 100
        d0_raw ~ Gamma(10, 1.0);  d0 = d0_raw / 100

        T = eltype(c)  # to get dual data type for AD
        att   = Vector{T}(undef, N)
        denom = zeros(T, nfrom)
        preds = Vector{T}(undef, N)

        @inbounds for i in 1:N
            att[i] = tp[i] * (l + (1-l)/((di[i] + d0)^c))
            if normalize
                denom[fromdist[i]] += att[i]
                denom[fromdist[i]] += exp(b) * (fpt[fromdist[i]] * (l + (1 - l) /
                    ((radius[fromdist[i]] + d0)^c)))
            end
        end

        @inbounds for i in 1:N
            if normalize
               preds[i] = fp[i] * exp(a) * (att[i] / denom[fromdist[i]])
            else
                preds[i] = fp[i] * exp(a) * (att[i])
            end
        end

        flows .~ Poisson.(preds)
        return preds
    end

    mdl = model(flows, fromdist, fp, tp, di, nfrom, N,
                fpt, radius, normalize)
    lb = [-20.0, -100.0, 10.0, 0.0, 1.0]
    ub = [20.0, 10.0, 100.0, 99.0, 100.0]
    return (; mdl, lb, ub)
end
