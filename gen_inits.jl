gen_random_inits = function(Nages, ncoefs)
    inits = [
        rand(Normal(0.0, 1.0), Nages); #a prior Exp = 0
        rand(Uniform(1.5, 2.5), Nages); #c prior Exp = 2.5
        rand(Uniform(0.3, 2.0), Nages); #d0 prior Exp = 2.5
        rand(Uniform(0.25, 1), Nages); #dscale prior Exp = 1
        [1.5];                     # netterr prior Exp = 7.5
        [8.0];                    # logconst prior Exp = 0
        fill(0.0, Nages); #kd                prior Exp = 0
        fill(0.0, Nages * ncoefs) # desirecoefs prior Exp = 0
        ]
  return inits
end

gen_bounds = function(Nages, ncoefs)
    a_lb = - 5.5
    a_ub = 5.5
    c_lb = .05
    c_ub = 30.0
    d0_lb = 0.0
    d0_ub = 60.0
    dscale_lb = .02
    dscale_ub = 3.0
    ne_lb = .5
    ne_ub = 50.0
    lc_lb = -30.0
    lc_ub = 30.0
    kd_lb = -2.5
    kd_ub = 2.5
    cheby_lb = -10.0
    cheby_ub = 10.0
        
    lower = [fill(a_lb, Nages); 
            fill(c_lb, Nages); 
            fill(d0_lb, Nages); 
            fill(dscale_lb, Nages); # dscale in fractions of meddist
            [ne_lb, lc_lb];
            fill(kd_lb, Nages); 
            cheby_lb * ones(ncoefs * Nages)]

    upper = [fill(a_ub, Nages); 
            fill(c_ub, Nages); 
            fill(d0_ub, Nages); 
            fill(dscale_ub, Nages); #dscale in fractions of meddist
            [ne_ub, lc_ub];
            fill(kd_ub, Nages); 
            cheby_ub * ones(ncoefs * Nages)]
    return(lower, upper)
end


show_inits = function(Nages, ncoefs) 
    names  = [fill("a", Nages); fill("c", Nages); fill("d0", Nages); fill("dscale", Nages); 
              "netterr"; "logconst"; fill("kd", Nages); fill("desire", Nages * ncoefs)]

    inits = gen_random_inits(Nages, ncoefs)
    lb = gen_bounds(Nages, ncoefs, -10, 10)[1]
    up = gen_bounds(Nages, ncoefs, -10, 10)[2]
    dt = DataFrame(names = names, lower = lb, inits = inits, upper = up)
    return dt
end
