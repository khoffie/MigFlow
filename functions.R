read_data <- function(mod_name) {
    p_preds <- "predictions"
    p_coefs <- "fitted_models"

    gen_filenames <- function(mod_name) {
    f_preds <- paste0("FlowDataPreds", mod_name, ".csv")
    f_coefs <- paste0("opti", mod_name, ".csv")
    nms <- list(f_preds = f_preds, f_coefs = f_coefs)

    f2 <- file.exists(file.path(p_coefs, nms[[2]]))
    if(file.exists(file.path(p_preds, nms[[1]])) == FALSE) {
        stop(sprintf("%s does not exist", nms[[1]]))
    }
    if(file.exists(file.path(p_coefs, nms[[2]])) == FALSE) {
        stop(sprintf("%s does not exist", nms[[2]]))
    }
    return(nms)
    }

    check_creation_time <- function(mod_name) {
        f1 <- file.path(p_preds, nms[[1]])
        f2 <- file.path(p_coefs, nms[[2]])
        message(sprintf("File %s created %s ",nms[1], file.info(f1)$ctime))
        message(sprintf("File %s created %s ",nms[2], file.info(f2)$ctime))
    }

    nms <- gen_filenames(mod_name)
    check_creation_time(mod_name)
    dt <- data.table::fread(file.path(p_preds, gen_filenames(mod_name)[1]))
    dt[, model := mod_name]
    coefs <- data.table::fread(file.path(p_coefs, gen_filenames(mod_name)[2]))
    coefs[, model := mod_name]
    out <- list(preds = dt, coefs = coefs)
    return(out)
}

calculate_net <- function(dt, col, by, o = "fromdist", d = "todist",
                           long = FALSE, type = NULL) {
    dt_in <- dt[, .(inflow = sum(get(col))), keyby = c(d, by)]
    dt_out <- dt[, .(outflow = sum(get(col))), keyby = c(o, by)]
    dt_in[dt_out, outflow := i.outflow, on = c(setNames(o, d), by)]
    dt_in[, net := inflow - outflow]
    setnames(dt_in, d, "region")
    if(!is.null(type)) {
        dt_in[, type := type]
        id_vars <- c("type", "region", by)
    }
    if(is.null(type)) {
        id_vars <- c("region", by)
    }
    if(long == TRUE) {
        dt_in <- melt(dt_in, id.vars = id_vars)
    }
    return(dt_in)
}


order_models <- function(dt) {
    model_names <- dt[, unique(model)]
    lvls <- data.table(models = model_names)
    lvls[, id := tstrsplit(model_names, "_")[[2]]]
    lvls <- lvls[order(as.numeric(id)), .(models)]
    dt[, model := factor(model, levels = lvls[, models])]
    return(NULL)
}

rec_ages <- function(dt) {
  agegroup <- NULL
  lbls <- data.table(old = c("unter18", "über65"),
                     new = c("below18", "above65"))
  dt[lbls, agegroup := i.new, on = .(agegroup = old)]
  return(NULL)
}

add_mising_flows <- function(flows, regions, agegroups, years) {
    ### in flows data all 0 flows are missing. We add them now to make
    ### sure all origins have the same destinations for all age groups and
    ### vice versa
    all_keys <- CJ(fromdist = regions,
                   todist = regions,
                   agegroup = agegroups,
                   year = years)
    setkeyv(all_keys, colnames(all_keys))
    setkeyv(flows, colnames(all_keys))
    flows <- flows[all_keys]
    flows[is.na(flows), flows := 0]
    return(flows)
}

fit_gravity <- function(dt, offset) {
    if(offset == TRUE) {
        fit <- glm(actual ~ offset(log(frompop)) +
                       offset(log(topop)) + log(distance),
                   family = poisson, data = dt)
    }
    if(offset == FALSE) {
        fit <- glm(actual ~ log(frompop) +
                       log(topop) + log(distance),
                   family = poisson, data = dt)
    }
    return(fit)
}

fit_zeroinfl <- function(dt) {
    fit <- pscl::zeroinfl(actual ~ log(frompop) + log(topop) + log(distance),
                          dist = "poisson", link = "logit", data = dt #, start = c(1, 1, -2)
                          )
    return(fit)
}

logistic <- function(x) {
    y <- 1 / (1 + exp(-x))
    return(y)
}

