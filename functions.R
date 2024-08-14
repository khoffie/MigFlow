read_data <- function(path, file) {
    if(grepl("prediction", path)) {
        mod <- gsub("FlowData", "", file)
        mod <- gsub("Preds", "model", mod)
    }
    if(grepl("models", path)) {
        mod <- gsub("opti_", "", file)
    }
    mod <- gsub(".csv", "", mod)
    dt <- data.table::fread(file.path(path, file))
    dt[, model := mod]
    return(dt)
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

plot_fit <- function(dt, x, y, th) {
    plt <- ggplot(dt[preds > th], aes({{x}}, {{y}})) +
        geom_hline(yintercept = 0) +
        geom_point(pch = ".") +
        geom_smooth(se = FALSE) +
        facet_wrap(vars(model, agegroup), scale = "free") +
        ggtitle(sprintf("Individual flows for preds > %s", th)) +
        theme_minimal()
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
