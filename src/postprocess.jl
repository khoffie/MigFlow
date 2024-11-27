function postprocess(; first = 50, path = nothing, render_plots = true, render_doc = true,
                     denstype = "all")
    file = "./writeup/juliaout_path.txt"
    if !isnothing(path); write(file, path); end
    if isnothing(path); path = readline(file); end
    mkpath(joinpath(path, "plots"))

    if render_plots
        saveparams(path, "germchain", kd)
        saveparams(path, "germchain", gravity)
        saveparams(path, "germchain", desire)
        saveparams(path, "germchain", densitychains, type = denstype)
        println("Plots generated")
    else
        println("Plots not generated")
    end
    ### confusingly compilereport uses path as in writeup/juliaout_path.txt
    compilereport(render_doc)
end

function compilereport(render_doc)
    # report.Rmd reads julia_output_path from file = "./writeup/juliaout_path.txt"
    if render_doc
        f = "./writeup/_main.Rmd"
        if isfile(f);  rm(f); end
        R"helpeR::render_doc('./writeup', 'report.Rmd')"
        print("Report generated")
    else
        println("Report not generated")
    end
end

function saveparams(path, pattern_in, fun; kwargs...)
    fin, fout = fileinout(path, pattern_in, fun)
    for i in eachindex(fin)
        chain = deserialize(joinpath(path, fin[i]))
        p = fun(chain; kwargs...)
        savefig(joinpath(path, "plots", "$(fout[i]).pdf"))
    end
end

function densitychains(chain; type = "best", densmin = 0, densmax = 5000, title = nothing)
    function plotdensity(chain, densmin, densmax)
        pars = OrderedDict(zip(string.(chain.value.axes[2]), chain.value.data[end, :, 1]))
        kds = [k for k in keys(pars) if contains(k, "kd")]
        kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),
                    getindex.(Ref(pars), kds) ./ 10)
        p = Plots.heatmap(kdfun, colorbar = false, ticks = false)
        return p
    end
    if type == "all"
        plts = [plotdensity(chain[:, :, i], densmin, densmax) for i in 1 : size(chain)[3]]
        p = Plots.plot(plts..., plot_title = title == nothing ? "" : title)
    elseif type == "best"
        bestchain = chain[:, :, findmax(chain[end, :lp, :])[2]]
        p = plotdensity(bestchain, densmin, densmax)
    end
    return p
end

function gravity(chain)
    ps = [:lp, :a, :c, :d0, :dscale, :e, :ktopop]
    p = moreparams(chain, ps, first)
    return p
end

function kd(chain)
    ps = [Symbol("kd[$i]") for i in [1, 5, 10, 15, 20, 30]]
    p = moreparams(chain, ps, first)
    return p
end

function desire(chain)
    ps = [Symbol("desirecoefs[$i]") for i in [1, 5, 10, 15, 20, 30]]
    p = moreparams(chain, ps, first)
    return p
end

function moreparams(chain, ps, first)
    out = Dict{Symbol, Any}()
    for p in ps
        out[p] = param(chain, p, first)
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

# function param(chain, symbol, first = 50)
#     x = range(first, size(chain)[1], step = 1)
#     p = Plots.plot(x, chain[symbol].data[first : end, :],
#                    xlab = string(symbol), label = "")
#     return p
# end

function param(chain, symbol, first)
    p = Plots.plot(chain[symbol].data[50 : end, :],
                   xlab = string(symbol), label = "")
    return p
end
