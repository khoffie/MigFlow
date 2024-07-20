#### My main issue is that in run_optimx() getting the predictions
#### often does not work. Then I need to fetch the optimal values and
#### call exp_moves() later. Also for different models I need
#### different exp_moves() but then I need to tell objective() to call
#### another exp_moves(). I would like to use the same code for
#### different models and simply loop over age groups in the data to
#### get predictins and coefs.


run_optimx <- function(dt, params, exp_moves) {
    objective = function(params,actual,...) {
        preds <- exp_moves(params, ...)
        out <- -sum(dpois(actual, preds, log = TRUE))
    }
    out <- optimx::optimx(par = params,
                          fn = objective,
                          actual = dt[, actual],
                          pop_o = dt[, pop_o],
                          pop_d = dt[, pop_d],
                          dist = dt[, dist],
                          method = "nlminb")
    optims <- ret_pars(out, pars = c("a", "b", "c", "d0", "e"))
    preds <- data.table(preds = exp_moves(params = optims,
                                          pop_o = dt[, pop_o],
                                          pop_d = dt[, pop_d],
                                          dist = dt[, dist]))
    out <- list(results = out, predictions = preds)
    return(out)
}

objective <- function(params, actual, ...) {
    preds <- exp_moves(params, ...)
    out <- -sum(dpois(actual, preds, log = TRUE))
  return(out)
}

exp_moves <- function(params, pop_o, pop_d, dist) {
    p <- params
    a <- p[["a"]]
    b <- p[["b"]]
    c <- p[["c"]]
##    d <- p[["d"]]
    d0 <- p[["d0"]]
    e <- p[["e"]]
##    preds <-  pop_o*pop_d*((a/1e4)*(1+b/(dist/(100 * e) + d0)^c))
    preds <-  pop_o*pop_d*((a/1e4)*(1+b/(dist/100 + d0)^c))
    return(preds)
}

ret_pars <- function(out, pars = c("a", "b", "c", "d", "d0"), digits = NULL) {
    p <- out$results[pars]
    if(!is.null(digits)) {
        p <- round(p, digits)
    }
    return(p)
}

##################### Workflow ###################

inis <- c(a = 1, b = 1, c = 1, d = 1, d0 = .05, e = 1)
out <- run_optimx(real, inis, objective)
optims <- ret_pars(out, c("a", "b", "c", "d0", "e"), NULL)
real[, preds_d0 := exp_moves(optims, pop_o, pop_d, dist)]

########################## not important ########################################

## ## ps <- make_paths(project_path = "~/Diss")
## ## real <- fread(file.path(ps$clean_flows_di, "districts_2000_2017_ger.csv"))
## ## dist <- fread(file.path(ps$dist, "distances_di.csv")) ## This
## ## real[dist, dist := i.distance_pos, on = .(origin, destination)]

## real <- ana$flows[year == 2017 & age_group != "all",
##                   .(origin, destination, age_group,
##                     dist = distance_pos, actual = od_flow)]
## real[age_dm[year == 2017], pop_o := i.pop, on = .(origin = region, age_group)]
## real[age_dm[year == 2017], pop_d := i.pop, on = .(destination = region, age_group)]
## real[age_group == "über65", age_group := "above65"]
## real[age_group == "unter18", age_group := "below18"]
## setnames(real, c("age_group", "origin", "destination"),
##          c("agegroup", "fromdist", "todist"))
## setcolorder(real, c("fromdist", "todist", "agegroup", "dist", "actual"))

## real[, pop_o := pop_o / 1e3]
## real[, pop_d := pop_d / 1e3]

## dm <- create_design_mat(aux$ink_di, 365, "Kreise", 2017)
## real[dm, bip_o := i.X365, on = .(fromdist = AGS)]
## real[dm, bip_d := i.X365, on = .(todist = AGS)]
## real[, gdp_ratio := bip_d / bip_o]

## ages <- c("below18", "18-25", "30-50")

## real <- real[agegroup == "25-30"]





### A Turing.jl sketch of a model

using Turing

@model migration1(flows,fromdist,todist,frompop,topop,distance,gdpc,agegroup,Ndist,meddist)

    a ~ Gamma(5.0,1.0/4.0)
    b ~ Gamma(5.0,1.0/4.0)
    c ~ Gamma(5.0,1.0/4.0)
    d0 ~ Gamma(5.0,1.0/4.0)

    preds = [frompop * topop * a/1000.0 * (1.0 + b/(dist/meddist + d0)^c) for _ in flows]
    flows ~ arraydist([Poisson(p) for p in preds])
    
end

function desir(d1,d2,d3,d4,age,from,to)
    if age == 1
        d1[to]/d1[from]
    elseif age == 2
        d2[to]/d2[from]
    elseif age == 3
        d3[to]/d3[from]
    elseif age == 4
        d4[to]/d4[from]
    end
end

@model migration2(flows,fromdist,todist,frompop,topop,distance,gdpc,agegroup,Ndist,meddist)

    a ~ Gamma(5.0,1.0/4.0)
    b ~ Gamma(5.0,1.0/4.0)
    c ~ Gamma(5.0,1.0/4.0)
    d0 ~ Gamma(5.0,1.0/4.0)
    desir1 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir2 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir3 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desir4 ~ arraydist([Gamma(5.0,1.0/4.0) for i in 1:Ndist])
    desires = [desir(desir1,desir2,desir3,desir4,agegroup[i],fromdist[i],todist[i]) for i in 1:length(flows)]
    preds = [frompop * topop * a/1000.0 * (1.0 + b/(dist/meddist + d0)^c) * desires ] 
    flows ~ arraydist([Poisson(p) for p in preds])
end

