using Plots, Distributions, StatsPlots

S = 10^5
N = 10^4

round2(x) = round(x, digits = 2)

function makeplots(p, S)
    N = length(p)
    xpb = rand(PoissonBinomial(p), S);
    xpois = rand(Poisson(N * mean(p)), S);
    xbinom = rand(Binomial(N, mean(p)), S);
    tvpois = round(.5sum(abs.(xpb .- xpois)), digits = 2)
    tvbinom = round(.5sum(abs.(xpb .- xbinom)), digits = 2)
    if sum(p) < 30
        f = Plots.histogram
        f! = Plots.histogram!
        al = .4
    else
        f = Plots.density
        f! = Plots.density!
        al = 1.0
    end
    p = f(xpb, label = "PB", alpha = al)
    p = f!(xpois, label = "Poisson", alpha = al)
    p = f!(xbinom, label = "Binomial", alpha = al)
    display(plot(p, title = "Error Pois $(tvpois) Error Binom $(tvbinom)"))
    return xpb, xpois, xbinom
end

p = vcat(rand(Uniform(0, 0.1), Int(N / 2)), rand(Uniform(0.0, 1.0), Int(N / 2)))
p = vcat(rand(Uniform(0.9, 0.95), Int(N / 2)), rand(Uniform(0.8, 1.0), Int(N / 10)))
p = vcat(rand(Uniform(0.001, 0.005), Int(N / 2)), rand(Uniform(0.9, 1.0), Int(N / 10)))

xpb, xpois, xbinom = makeplots(p, S);

x = 1:N

pdfpb = pdf.(PoissonBinomial(p), x)
pdfb = pdf.(Binomial(length(p), mean(p)), x)
pdfp = pdf.(Poisson(sum(p)), x)

plot(pdfpb, xlim = (2500, 3000))
plot!(pdfb)
plot!(pdfp)

tverror(p,q ) = round2(.5sum(abs.(p .- q)))

tverror(pdfpb, pdfb)
tverror(pdfpb, pdfp)
