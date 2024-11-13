
function param(chain, symbol)
    p = Plots.plot(chain[symbol], xlab = string(symbol), label = "")
    return p
end

function moreparams(chain, ps)
    out = Dict{Symbol, Any}()
    for p in ps
        out[p] = param(chain, p)
    end
    p = Plots.plot([out[i] for i in ps]...)
    return p
end

function gravity(chain)
    ps = [:a, :c, :d0, :dscale, :e, :ktopop]
    p = moreparams(chain, ps)
    return p
end

function kd(chain)
    ps = [Symbol("kd[$i]") for i in [1, 5, 10, 15, 20, 30]]
    p = moreparams(chain, ps)
    return p
end

function desire(chain)
    ps = [Symbol("desirecoefs[$i]") for i in [1, 5, 10, 15, 20, 30]]
    p = moreparams(chain, ps)
    return p
end

fns = readdir(path)
occursin.("chain", fns)

function allparams(path)
    
