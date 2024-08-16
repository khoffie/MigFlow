 function gen_random_inits(Nages, ncoefs)
     inits = [
         rand(Normal(0.0, 1.0), Nages); #a prior Exp = 0
         rand(Uniform(1.5, 2.5), Nages); #c prior Exp = 2.5
         rand(Uniform(1.0, 10.0), Nages); #d0 prior Exp = 2.5
         rand(Uniform(0.25, 1), Nages); #dscale prior Exp = 1
         [1.5];                     # netterr prior Exp = 7.5
         [5.0];                    # logconst prior Exp = 0
         fill(0.0, Nages); #kd                prior Exp = 0
         fill(0.0, Nages * ncoefs) # desirecoefs prior Exp = 0
         ]
   return inits
 end

function gen_fixed_inits(Nages, ncoefs)   
## handpicked inits from manualinits.ipynb
    #= as = [.7, 2.3, 2.5, .1, -.3, -.3]
   ##as = [-10, -10, -10, -10, -10, -.3] 
    cs = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
    d0s = [.1, .1, .1, .1, .1, .1]
    dscales = [.3, .3, .5, .5, .3, .3]
    logconst = -5
   ##as = fill(-2, Nages) =#

   inits = [
    [.7, 2.3, 2.5, .1, -.3, -.3];
    [2.5, 2.5, 2.5, 2.5, 2.5, 2.5];
    [.1, .1, .1, .1, .1, .1];
    [.3, .3, .5, .5, .3, .3];
   [1.5];
   [5.0];
   fill(0.0, Nages);
   fill(0.0, Nages * ncoefs)
   ]
   return inits
end


function gen_bounds(Nages, ncoefs)
    a_lb = - 1.0
    a_ub = 3.0
    c_lb = 1.5
    c_ub = 3.5
    d0_lb = 0.0
    d0_ub = 1.0
    dscale_lb = .1
    dscale_ub = .5
    ne_lb = .5
    ne_ub = 10.0
    lc_lb = 4.0
    lc_ub = 6.0
    kd_lb = -0.5
    kd_ub = 0.5
    

    cheby_lb = -0.01
    cheby_ub = 0.01
        
    lower = [
        ##fill(a_lb, Nages); 
        [0.0, 2.0, 2.0, -.5, -1.0, -1.0]
        fill(c_lb, Nages); 
        fill(d0_lb, Nages); 
        fill(dscale_lb, Nages); # dscale in fractions of meddist
        [ne_lb, lc_lb];
        fill(kd_lb, Nages); 
        cheby_lb * ones(ncoefs * Nages)]

    upper = [
    ##    fill(a_ub, Nages); 
        [1.0, 3.0, 3.0, 1.0, 1.0, 1.0];
        fill(c_ub, Nages); 
        fill(d0_ub, Nages); 
        fill(dscale_ub, Nages); #dscale in fractions of meddist
        [ne_ub, lc_ub];
        fill(kd_ub, Nages); 
        cheby_ub * ones(ncoefs * Nages)]
    return(lower, upper)
end

 function gen_inits_bounds(; Nages, ncoefs, type, opts_f, show) 
    names  = [fill("a", Nages); fill("c", Nages); fill("d0", Nages); fill("dscale", Nages); 
              "netterr"; "logconst"; fill("kd", Nages); fill("desire", Nages * ncoefs)]
    if type == "random"
        inits = gen_random_inits(Nages, ncoefs)
    elseif type == "fixed"
        inits = gen_fixed_inits(Nages, ncoefs)
    elseif type == "opts"
        inits = CSV.read(opts_f, DataFrame)[:, "values"]
    end
    lb = gen_bounds(Nages, ncoefs)[1]
    up = gen_bounds(Nages, ncoefs)[2]
    ib = DataFrame(names = names, lowers = lb, inits = inits, uppers = up)
    if show == true
        ## better print one of each param
        @printf("Inits and bounds %s\n", ib[1:40, :])
    end
    return ib
end
