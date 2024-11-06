library(data.table)

correct_flows <- function(f, c) {
    ## if speed is an issue. Following three lines should be done outside the function
    fm <- as.matrix(dcast(f, origin ~ destination, value.var = "flow", fill = 0)[, -1])
    cm <- as.matrix(dcast(c, ags_old ~ ags_new, value.var = "conv_p", fill = 0)[, -1])
    rownames(cm) <- c[, unique(ags_old)]

    fnewm <- t(cm) %*% fm %*% cm
    
    fnew <- melt(data.table(fnewm, keep.rownames = TRUE), id.vars = "rn")
    setnames(fnew, c("rn", "variable", "value"), c("origin", "destination", "flow"))
    fnew[, flow := as.integer(round(flow))]
    return(fnew)
}

## according correct.csv 3201 and 3253 did not exist in 2001
## anymore. In flows they only appear as origin, not as
## destination. Probably they are coded wrongly. For now, I simply
## remove them.
dt_flows <- fread("./data/raw/flows_districts_2000_2017_ger.csv")
correct <- data.table::fread(file.path("~/Diss/inst/extdata/clean/corrections", "districts_19.csv"))

dt_flows <- dt_flows[! dt_flows[year == 2001 & origin %in% c(3201, 3253)],
                     on = .(origin, destination, year, age_group)]
dt_new <- dt_flows[, correct_flows(.SD, correct[year == .BY$year]), keyby = .(year, age_group)]
