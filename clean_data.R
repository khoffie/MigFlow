library(data.table)
library(helpeR)
library(sf)

shp <- data.table(sf::read_sf("./data/raw/shapes/districts_ext.shp"))

clean_pop <- function(path, file = c("12411-03-02-4.csv", "12411-03-03-4.csv"), shp) {
    file <- match.arg(file)
    dt <- read_pop(path, file)
    new_cols <- c("year", "region", "region_type", "age_group",
                  "all_b", "all_m", "all_f",
                  "ger_b", "ger_m", "ger_f",
                  "for_b", "for_m", "for_f")
    setnames(dt, colnames(dt), new_cols)

    dt <- rec_ages(dt, file)

    dt <- dt[-c(1,2), ] ## colnames in rows
    dt <- dt[region_type != "Deutschland"] ## rm since interested in districts
    dt[, year := as.integer(substring(year, 7, 10))]
    dt[, region := as.integer(region)]
    dt[region == 2, region := 2000] ## Hamburg
    dt[region == 11, region := 11000] ## Berlin

    nc <- new_cols[-c(1:4)]
    dt[, c(nc) := lapply(.SD, as.integer), .SDcols = nc]
    
    dt <- dt[, .(year, region, age_group, all_b, ger_b, for_b)]
    dt <- dt[region %in% shp[, AGS]] # filter all non districts

    dt <- dt[, c("all_b", "ger_b", "for_b") := lapply(.SD, sum),
             keyby = .(year, region, age_group), .SDcols = c("all_b", "ger_b", "for_b")]

    dt <- dt[, .SD[1], keyby = .(year, region, age_group)]
    dt <- dt[, .(region, year, age_group, all = all_b, german = ger_b, foreign = for_b)]

    return(dt)
}

rec_ages <- function(dt, file) {
    last <- ifelse(file == "12411-03-02-4.csv", 2, 6)
    rec <- dt[, .(old_age = unique(age_group))]
    a <- c("unter18", "18-25", "25-30", "30-50", "50-65", "über65" )
    new_ages <- c(NA, rep(a[1], 5), rep(a[2], 2), a[3], rep(a[4], 4), rep(a[5], 3), rep(a[6], last), "all")
    rec[, new_age := new_ages]
    dt[rec, age_group := i.new_age, on = .(age_group = old_age)]    
    return(dt)
}

read_pop <- function(path, file) {
    read_pop_ <- function(path, file) {
        dt <- fread(file.path(path, file), skip = 6,
                    encoding = "Latin-1", na.strings = c("-", "", "."),
                    strip.white = TRUE, fill = FALSE)
        return(dt)
    }
    ## first linenumbers not part of data
    cond <- ifelse(file == "12411-03-02-4.csv", "184005", "153877")
    cond <- sprintf("Stopped early on line %s", cond)
    dt <- withCallingHandlers({
        read_pop_(path, file)
    }, warning = function(w) {
        ## data ends there and only meta information follows
        if (grepl(cond, conditionMessage(w))) {
            invokeRestart("muffleWarning")
        } else {
            warning(w)
        }
    })
    return(dt)
}


dt1 <- clean_pop("./data/raw", "12411-03-02-4.csv", shp)
dt2 <- clean_pop("./data/raw", "12411-03-03-4.csv", shp)

dt <- rbind(dt1, dt2)
mis <- dt[year >= 2000 & age_group == "all", sum(is.na(german)), keyby = .(region)]
tail(mis[order(V1)], 20)

dt[region == 3159 & age_group == "all"]

shp[AGS == 3159]
