function postprocess(path = nothing)
    if isnothing(path); path = readline("./writeup/juliaout_path.txt"); end
    mkpath(joinpath(path, "plots"))
    saveparams(path, "germchain", kd)
    saveparams(path, "germchain", gravity)
    saveparams(path, "germchain", desire)
    saveparams(path, "germchain", densitychains)
end

function saveparams(path, pattern_in, fun)
    fin, fout = fileinout(path, pattern_in, fun)
    for i in eachindex(fin)
        chain = deserialize(joinpath(path, fin[i]))
        p = fun(chain)
        savefig(joinpath(path, "plots", "$(fout[i]).pdf"))
    end
end

function densitychains(chain; densmin = 0, densmax = 5000, title = nothing)
    function plotdensity(chain, densmin, densmax)
        pars = OrderedDict(zip(string.(chain.value.axes[2]), chain.value.data[end, :, 1]))
        kds = [k for k in keys(pars) if contains(k, "kd")]
        kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
                    getindex.(Ref(pars), kds) ./ 10)
        p = Plots.surface(kdfun, colorbar = false, ticks = false)
        return p
    end
    plts = [plotdensity(chain[:, :, i], densmin, densmax) for i in 1 : size(chain)[3]]
    p = Plots.plot(plts..., plot_title = title == nothing ? "" : title)
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

function moreparams(chain, ps)
    out = Dict{Symbol, Any}()
    for p in ps
        out[p] = param(chain, p)
    end
    p = Plots.plot([out[i] for i in ps]...)
    return p
end

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
