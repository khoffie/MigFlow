function fitmodel(path, alldata, sampler, n_samples, thinning, name)
    d = alldata
    chn = @time(Turing.sample(d.model, sampler, MCMCThreads(),
                                n_samples, 4, thinning = thinning,
                                initial_params = fill(d.vals.inits, 4),
                                verbose = true, progress = true))
    if contains(string(sampler), "Slice")
        ## lp values for slice are wrong
        chn[:, :lp, :] = logprob(d.model, chn)
    end
    Serialization.serialize(joinpath(path, "chains/$(name)"), chn)
    return chn
end

# function fitmodel(alldata, sam, outpath, thinning, show_plot = true)
#     alldata, optis, chn = @time(runsampling(alldata.model, alldata, sam,
#                                     alldata.vals.params,
#                                     fill(alldata.vals.inits, 4);
#                                     chainout = outpath,
#                                     nchains = 4,
#                                     nsamples = 100,
#                                     thinning = thinning,
#                                             paramtype = "last"))
#     if show_plot
#         display(plot(gravity(chn),
#                      plotfit(alldata.flows), layout = (2, 1)))
#     end

#     return alldata, chn
# end

maxlp(chain) = round(maximum(chain[end, :lp, :]))

function gravity(chain, first = 1, first_lp = 50, title = "")
    ps = [:lp, :a, :c, :d0, :dscale, :e, :ktopop]
    ps = intersect(ps, chain.name_map[1])
    ps = vcat(:lp, ps)
    plot([plotparam(chain, p, first, first_lp) for p in ps]...,
         plot_title = title)
end

function plotparam(chain, symbol, first, first_lp)
    if symbol == :lp first = first_lp end
    xrange = first : size(chain)[1]
    yrange = chain[symbol].data[first : end, :]
    p = Plots.plot(xrange, yrange, xlab = string(symbol), label = "")
    return p
end
