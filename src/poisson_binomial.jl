round2(x) = round(x, digits = 2)
tverror(p,q ) = round2(.5sum(abs.(p .- q)))

function approx_plot(ps, xlim = nothing)
    p = eval(Meta.parse(ps))
    x = 0:10length(p)
    pb = pdf.(PoissonBinomial(p), x)
    binom = pdf.(Binomial(length(p), mean(p)), x)
    pois = pdf.(Poisson(sum(p)), x)
    tvb = tverror(pb, binom)
    tvp = tverror(pb, pois)

    # if sum(p) < 30
    #     f = Plots.histogram
    #     f! = Plots.histogram!
    #     al = .4
    # else
        f = Plots.plot
        f! = Plots.plot!
        al = 1.0
##    end

    args = Dict(:label => "PB", :alpha => al)
    if xlim == nothing
        args[:xlim] = (.5sum(p), 2sum(p))
    else
        args[:xlim] = xlim
    end

    p = f(pb; args...)
    p = f!(pois, label = "Poisson", alpha = al)
    p = f!(binom, label = "Binomial", alpha = al)
    p = plot(p, title = "$(ps) \nError Pois $(tvp) Error Binom $(tvb)", titlefontsize = 8)
    display(p)
    return pb, pois, binom, p
end
