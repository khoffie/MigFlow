using Pkg
Pkg.activate(".")



using ApproxFun, StatsPlots, Optimization, OptimizationOptimJL, Distributions,
    CSV, DataFrames

function chebpoly(coefs,xrange)
    fun = Fun(Chebyshev(xrange),coefs)
end

plot(chebpoly([2.0,1.0,.2,.1],0..20); xlim=(0,20))


dcp = chebpoly([2.0,1.0,.2,.1,.03],0..20)
dcp2 = dcp*dcp
cp = Integral()*dcp2
cp = cp-cp(0)
plot(cp)


function expmoves(frompop,topop,dist,a,b,c,d0,gdprat,gdpfun)
    frompop * topop * (a/10000) * (1+b/(dist/100 + d0))^c * gdpfun(gdprat)
end

function makegdpfun(coefs)
    dgdpfunsqrt = Fun(Chebyshev(p.gdprange),coefs)
    dgdpfun = dgdpfunsqrt*dgdpfunsqrt;
    gdpfun = Integral()*(dgdpfun) #integrate a function squared... results in a non-decreasing function
    gdpfun = gdpfun - gdpfun(0.0) ## this is now a function that's 0 at 0 and non-decreasing

    return (gdpfun,dgdpfun)
end


function objfun(x,p)
    data = p.dat
    a = x[1]
    b = x[2]
    c = x[3]
    d0 = x[4]
    
    gdpfun,dgdpfun = makegdpfun(x[5:10])

    sum = 0.0
    for r in eachrow(data)
        gdprat = r.gdpto / r.gdpfrom # needs to be gdp/capita of from and to district
        sum = sum + logpdf(Poisson(expmoves(r.frompop,r.topop,r.dist,a,b,c,d0,gdprat,gdpfun)),
                                r.actual)

    end
    # penalized negative log likelihood (ie. incorporate prior info about the value and derivative at x=1.0)
    return -sum + 1000.0*(gdpfun(1.0)-1.0)^2 + 1000.0*(dgdpfun(1.0))^2
end


### to run this, we need input data, so you should 
### write out a data frame that has the following columns:
### from,to,actual,dist,frompop,topop,gdpfrom,gdpto
### from,to are the district ids, actual is actual migration, dist = distance in km, 
### frompop = population of from distric, topop = population of to district
### gdpfrom = gdp per CAPITA in from district, gdpto = gdp per CAPITA in to district

function run()

    thedf = CSV.read("data/optimizationdat.csv",DataFrame())


    gdprange = 0.0 .. ceil(maximum(thedf.gdpto ./ thedf.gdpfrom))

    optfun = OptimizationFunction(objfun,Optimization.AutoForwardDiff())

    prob = OptimizationProblem(optfun,[1.0,5.0,1.5,0.2,1.0,0.0,0.0,0.0,0.0,0.0],(dat=thedf,gdprange = gdprange))

    soln = solve(prob,BFGS())

    u = soln.u
    gdpfun = makegdpfun(u[5:10])
    a = u[1]
    b = u[2]
    c = u[3]
    d0 = u[4]
    ydat = [log(r.actual / expmoves(r.frompop,r.topop,r.dist,a,b,c,d0,) )]
    xdat = thedf.dist 

    p = plot(xdat,ydat)
    display(p)
    
end

run()

