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

make_map <- function(netm, age) {
    plt <- ggplot(set_geom(netm[agegroup == age], F)) +
        geom_sf(aes(fill = value)) +
        facet_wrap(~variable) +
        ggtitle(sprintf("Actual and predicted nmr (normalized by pop of Germans) for age group %s", age)) +
        theme_minimal() 
    return(plt)
}

make_net_plot <- function(net, scales = "free") {
    plt <- ggplot(net, aes(preds, actual)) +
        geom_hline(yintercept = 0, col = "blue", linewidth = .2) +
        geom_vline(xintercept = 0, col = "blue", linewidth = .2) +        
        geom_point(alpha = .5) +
        facet_wrap(vars(agegroup, year), scales = scales) +
        theme_minimal() +
        ggtitle("Net migration rate(400 regions, normalized by pop of Germans)")
    return(plt)
}

plot_fit <- function(dt, x, y, th_min, th_max = NULL) {
    main <- sprintf("Individual flows for preds > %s", th_min)
    plt <- ggplot(dt[preds > th_min], aes({{x}}, {{y}})) +
        geom_hline(yintercept = 0) +
        geom_point(pch = ".", alpha = .3) +
        geom_smooth(se = FALSE) +
        facet_wrap(vars(model, agegroup), scale = "free") +
        ggtitle(main) +
        theme_minimal()
    if(is.null(th_max) == FALSE) {
        main <- sprintf(paste(main, "and < %s"), th_max)
        plt <- plt + xlim(c(0, th_max)) +
            ggtitle(main)
    }
    return(plt)
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

make_net_map <- function(net) {
    plt <- ggplot(set_geom(net, F)) +
        geom_sf(aes(fill = value)) +
        facet_wrap(vars(model, type, agegroup)) +
        theme_minimal() 
    return(plt)
}
