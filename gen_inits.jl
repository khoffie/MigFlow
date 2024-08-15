gen_random_inits = function(Nages, ncoefs)
    inits = [
        rand(Normal(0.0, 1.0), Nages); #a
        rand(Uniform(1.5, 2.5), Nages); #c
        rand(Uniform(0.0, 0.03), Nages); #d0
        rand(Uniform(0.1, 1), Nages); #dscale
        [1.5, 8.0]; # neterr and logisticconst
        fill(0.0, Nages); #kd
        fill(0.0, Nages * ncoefs) # desirecoefs
        ]
  return inits
end

gen_bounds = function(Nages, ncoefs, cheby_lb, cheby_ub)
    a_lb = - 5.5
    a_ub = 5.5

    c_lb = .05
    c_ub = 10
    d0_lb = 0.0
    d0_ub = 60.0
    ne_lb = .5
    ne_ub = 50.0
    lc_lb = -30.0
    lc_ub = 30.0
    kd_lb = -2.5
    kd_ub = 2.5
        
    lower = [fill(a_lb, Nages); 
            fill(c_lb, Nages); 
            fill(d0_lb, Nages); 
            fill(0.02,Nages); # dscale in fractions of meddist
            [ne_lb, lc_lb];
            fill(kd_lb, Nages); 
            cheby_lb * ones(ncoefs * Nages)]

    upper = [fill(a_ub, Nages); 
            fill(c_ub, Nages); 
            fill(d0_ub, Nages); 
            fill(3.0,Nages); #dscale in fractions of meddist
            [ne_ub, lc_ub];
            fill(kd_ub, Nages); 
            cheby_ub * ones(ncoefs * Nages)]
    return(lower, upper)
end


check_inits = function(Nages, ncoefs) 
    names  = [fill("a", Nages); fill("c", Nages); fill("d0", Nages); fill("dscale", Nages); 
              "netterr"; "logconst"; fill("kd", Nages); fill("desire", Nages * ncoefs)]

    inits = gen_random_inits(Nages, ncoefs)
    lb = gen_bounds(Nages, ncoefs, -10, 10)[1]
    up = gen_bounds(Nages, ncoefs, -10, 10)[2]
    dt = DataFrame(names = names, lower = lb, inits = inits, upper = up)
    return dt
end
