function norm(data::NamedTuple; normalize = true, type)
    df = sort(data.df, :fromdist)
    districts = sort(data.districts, :distcode)
    ds = 100

    Y    = df.flows
    from = lc(df.fromdist)
    to = lc(df.todist)
    if type == "joint"
        A = df.frompop
    elseif type == "conditional"
        df2 = combine(groupby(data.df, :fromdist), :flows => sum)
        A = leftjoin(df, df2, on = :fromdist).flows_sum
    end
    P       = districts.pop ./ 153000 # median topop
    poporig = districts.pop
    D       = df.dist  ./ ds
    Ndist    = length(districts.distcode)
    N        = length(Y)
    radius   = fradius.(districts.pop, districts.density)
    data = (; Y, D, from, to, A,  P, poporig)

    @model function model(Y, from, to, A, P, D, Ndist, N,
                          radius, normalize)
        a      ~ Normal(-5, 1)
        b  ~ Gamma(1, 1);     ## b  = b_raw / 100
        c_raw  ~ Gamma(15, 0.2);  c  = c_raw / 10
        l_raw  ~ Gamma(10, 1.0);  l  = l_raw / 100
        d0_raw ~ Gamma(10, 1.0);  d0 = d0_raw / 100

        T = eltype(c)  # to get dual data type for AD
        att   = Vector{T}(undef, N)
        denom = zeros(T, Ndist)
        preds = Vector{T}(undef, N)

        @inbounds for i in 1:N
            att[i] = P[to[i]] * (l + (1-l)/((D[i] + d0)^c))
            if normalize
                denom[from[i]] += att[i]
                denom[from[i]] += exp(b) * (P[from[i]] * (l + (1 - l) /
                    ((radius[from[i]] + d0)^c)))
            end
        end

        @inbounds for i in 1:N
            if normalize
               preds[i] = A[i] * exp(a) * (att[i] / denom[from[i]])
            else
                preds[i] = A[i] * exp(a) * (att[i])
            end
        end

        Y .~ Poisson.(preds)
        return preds
    end

    mdl = model(Y, from, to, A, P, D, Ndist, N,
                radius, normalize)
    lb = [-20.0, -100.0, 10.0, 0.0, 1.0]
    ub = [20.0, 10.0, 100.0, 99.0, 100.0]
    return (; mdl, lb, ub, data)
end

fradius(P, ρ) = sqrt((P / ρ) / 2π)
lc(x) = levelcode.(categorical(x))
