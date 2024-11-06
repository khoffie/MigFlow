library(data.table)
library(MigStat)

correct <- data.table::fread(file.path("~/Diss/inst/extdata/clean/corrections", "districts_19.csv"))
correct[ags_old %in% c(3201, 3253)]


dt_flows <- fread("./data/raw/flows_districts_2000_2017_ger.csv")
correct_flows <- function(f, c) {
    fm <- as.matrix(dcast(f, origin ~ destination, value.var = "flow", fill = 0)[, -1])
    cm <- as.matrix(dcast(c, ags_old ~ ags_new, value.var = "conv_p", fill = 0)[, -1])
    rownames(cm) <- c[, unique(ags_old)]

    fnewm <- t(cm) %*% fm %*% cm
    
    fnew <- melt(data.table(fnewm, keep.rownames = TRUE), id.vars = "rn")
    setnames(fnew, c("rn", "variable", "value"), c("origin", "destination", "flow"))
    fnew[, flow := as.integer(round(flow))]
    return(fnew)
}

dt_new <- dt_flows[year > 2001, correct_flows(.SD, correct[year == .BY$year]), keyby = .(year, age_group)]
dt_new[, origin := as.integer(origin)]
dt_new[, destination := as.integer(as.character(destination))]

## Old approach
flows_new <- MigStat::correct_flows(dt_flows, correct)
dt_new[flows_new, old := i.flow, on = .(year, age_group, origin, destination)]
dt_new[origin != destination][flow != old]


f <- dt_flows[year == 2001 & age_group == "all"]
c <- correct[year == 2001]

setdiff(f[, unique(origin)], c[, unique(ags_old)])
f[origin == 3201]






fnew[origin != destination]
flows_new[year == 2000 & age_group == "18-25"][order(destination)]


flows_new <- correct_flows(dt_flows, correct[conv_p != 0])
dcast(flows_new, origin ~ destination)
dcast(flows, origin ~ destination, value.var = "flow")
MigStat:::correct_flow(flows, changes, ags_col = "destination")


############################ example #################################
f <- c(0, 100, 20, 20, 0, 100, 50, 10, 0)
fm <- matrix(c(0, 20, 50, 100, 0, 10, 20, 100, 0), ncol = 3)
o <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
d <- c("A", "B", "C", "A", "B", "C", "A", "B", "C")
flows <- data.table(origin = o, destination = d, flow = f, year = 2017)

c <- c(1, 0, 0, 1, 0, 0, .2, 0, .8)
cm <- matrix(c(1, 1, .2, 0, 0, 0, 0, 0, .8), ncol = 3)
t(cm) %*% fm %*% cm
fm

fm2 <- c(0, 10, 10, 100, 0,
         100, 0, 20, 10, 0,
         20, 100, 0, 10, 0,
         10, 20, 100, 0, 0,
         0, 0, 0, 0, 0)
fm2 <- matrix(fm2, ncol = 5, byrow = TRUE)
cm2 <- c(1, 0, 0, 0, 0, .2, .8, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
cm2 <- matrix(cm2, ncol = 5, byrow = TRUE)

new <- t(cm2) %*% fm2 %*% cm2
