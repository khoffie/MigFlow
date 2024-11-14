function fileinout(path, pattern_in, fun)
    funname = string(fun)
    function newname(files, new)
        fnnew = replace.(files, pattern_in => funname)
        fnnew = [f[begin : (end - 4)] for f in fnnew]
        return fnnew
    end
    fin = readdir(path)
    fin = fin[occursin.(pattern_in, fin)]
    fout = newname(fin, funname)
    return fin, fout
end

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
    ps = [:lp, :a, :c, :d0, :dscale, :e, :ktopop]
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

function saveparams(path, pattern_in, fun)
    fin, fout = fileinout(path, pattern_in, fun)
    for i in eachindex(fin)
        chain = deserialize(joinpath(path, fin[i]))
        p = fun(chain)
        savefig(joinpath(path, "plots", "$(fout[i]).pdf"))
    end
end

function allparams(path)
    saveparams(path, "germchain", kd)
    saveparams(path, "germchain", gravity)
    saveparams(path, "germchain", desire)
end

path = "manuscript_input/1000slice"
allparams(path)


fileinout(path, "densfun", density)

"germ" * "dens"
