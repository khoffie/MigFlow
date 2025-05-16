function norm(data::NamedTuple; normalize = true)
    flows    = data.flows
    from     = data.fromdist
    fp       = data.frompop
    tp       = data.topop ./ 153000 # median topop
    di       = data.dist  ./ data.distscale
    nfrom    = length(unique(from))
    N        = length(flows)
    fpt      = data.fpt
    radius   = data.radius
    fromfull = data.fromdistfull
    tpfull   = data.topopfull
    Nfull    = length(data.fromdistfull)
    difull   = data.distfull
    @model function model(flows, from, fp, tp, di, nfrom, N,
                          fpt, radius, normalize, fromfull, tpfull,
                          Nfull)
        # Priors
        a      ~ Normal(-5, 1)
        b      ~ Gamma(1, 1);
        c_raw  ~ Gamma(15, 0.2);  c  = c_raw / 10
        l_raw  ~ Gamma(10, 1.0);  l  = l_raw / 100
        d0_raw ~ Gamma(10, 1.0);  d0 = d0_raw / 100

        T = eltype(c)  # to get dual data type for AD
        att   = Vector{T}(undef, N)
        denom = zeros(T, nfrom)
        preds = Vector{T}(undef, N)

        @inbounds for i in 1:N
            att[i] = desirability(tp[i], di[i], l, d0, c)
        end

        # if normalize
        #     @inbounds for i in 1:Nfull
        #         denom[fromfull[i]] += desirability(tpfull[i], difull[i], l, d0, c)
        #         denom[fromfull[i]] += exp(b) *
        #             desirability(fpt[fromfull[i]], radius[fromfull[i]], l, d0, c)
        #     end
        # end

        if normalize
            @inbounds for i in 1:N
                denom[from[i]] += desirability(tp[i], di[i], l, d0, c)
                denom[from[i]] += exp(b) *
                    desirability(fpt[from[i]], radius[from[i]], l, d0, c)
            end
        end

        @inbounds for i in 1:N
            if normalize
               preds[i] = fp[i] * exp(a) * (att[i] / denom[from[i]])
            else
                preds[i] = fp[i] * exp(a) * (att[i])
            end
        end

        flows .~ Poisson.(preds)
        return preds
    end

    mdl = model(flows, from, fp, tp, di, nfrom, N,
                fpt, radius, normalize, fromfull, tpfull, Nfull)
    lb = [-20.0, -100.0, 10.0, 0.0, 1.0]
    ub = [20.0, 10.0, 100.0, 99.0, 100.0]
    return (; mdl, lb, ub)
end

desirability(P, D, ϕ, δ, γ) = P * (ϕ + (1 - ϕ) / (D + δ)^γ)
