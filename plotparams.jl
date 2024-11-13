
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

function lp(chain)
    p = moreparams(chain, [:lp])
    display(p)
    return p
end

function saveparams(path, fun)
    funname = string(fun)
    function newname(files, new)
        fnnew = replace.(files, "germchain" => funname)
        fnnew = [f[begin : (end - 4)] for f in fnnew]
        return fnnew
    end
    fns = readdir(path)
    fnchains = fns[occursin.("chain", fns)]
    fnnew = newname(fnchains, funname)
    for i in 1 : length(fnnew)
        chain = deserialize(joinpath(path, fnchains[i]))
        p = fun(chain)
        savefig(joinpath(path, "plots", "$(fnnew[i]).pdf"))
    end
end

function allparams(path)
    saveparams(path, kd)
    saveparams(path, gravity)
    saveparams(path, desire)
    saveparams(path, lp)
end
