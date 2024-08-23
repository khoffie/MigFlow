make_map <- function(netm, age) {
    main <- sprintf("Actual and predicted nmr (normalized by pop of Germans) for age group %s", age)
    plt <- ggplot(set_geom(netm[agegroup == age], F)) +
        geom_sf(aes(fill = value)) +
        facet_wrap(~variable) +
        ggtitle(main) +
        theme_minimal() 
    return(plt)
}

make_net_plot <- function(net, type = "bivariate", scales = "free") {
    normalized <- net[, sum(actual > 1) == 0] ## if TRUE, data is not normalized
    if(normalized == TRUE) {
        main <- "Net migration rate(400 regions, normalized by pop of Germans)"
    }
    if(normalized == FALSE) {
         main <- "Net migration"
    }
    if(type == "bivariate") {
        plt <- ggplot(net, aes(preds, actual)) +
        geom_hline(yintercept = 0, col = "blue", linewidth = .2) +
        geom_vline(xintercept = 0, col = "blue", linewidth = .2) +        
        geom_point(alpha = .5) 
    }
    if(type == "univariate") {
        plt <- ggplot(netm[type != "diff"], aes(value, fill = type)) +
            geom_density(alpha = .2) +
            facet_wrap(vars(agegroup), scales = "free") +
            ggtitle("Distribution of net migration rate") +
            theme_minimal()
    }
    plt <- plt +
        facet_wrap(vars(agegroup, year), scales = scales) +
        theme_minimal() +
        ggtitle(main)
    return(plt)
}

plot_fit <- function(dt, x, y, th_min, smooth = TRUE, th_max = NULL, p_sample = NULL) {
    main <- "Individual flows "
    if(!is.null(p_sample)) {
        dt <- dt[sample(1:.N, size = as.integer(p_sample * .N))]
        main <- sprintf("%s percent sample of individual flows", p_sample * 100)
    }
    main <- sprintf(paste(main, "for preds > %s"), th_min)
    plt <- ggplot(dt[preds > th_min], aes({{x}}, {{y}})) +
        geom_hline(yintercept = 0) +
        geom_point(pch = ".", alpha = .3) +
        facet_wrap(vars(agegroup, model), scale = "fixed") +
        ggtitle(main) +
        theme_minimal()
    if(is.null(th_max) == FALSE) {
        main <- sprintf(paste(main, "and < %s"), th_max)
        plt <- plt + xlim(c(0, th_max)) +
            ggtitle(main)
    }
    if(smooth == TRUE) {
        plt <- plt + geom_smooth(se = FALSE) 
    }
    return(plt)
}

make_net_map <- function(net) {
    plt <- ggplot(set_geom(net, F)) +
        geom_sf(aes(fill = value)) +
        facet_wrap(vars(model, type, agegroup)) +
        theme_minimal() 
    return(plt)
}

ggsavew <- function(w = 8, ...) {
    w <- w
    h <- w / 1.6
    ggsave(..., width = w, height = h)
}
