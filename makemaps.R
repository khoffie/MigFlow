library(data.table)
library(sf)
library(MigStat)
library(ggplot2)
library(patchwork)
source("functions.R")
path <- "./manuscript_input/bu"

read_output <- function(path, type) {
    read_agegroup <- function(path, file) {
        dt <- data.table::fread(file.path(path, file))
        dt[, agegroup := get_age(file)]
    return(dt)
}
    get_age <- function(file) {
        age <- data.table::tstrsplit(file, "_")[[2]]
        age <- data.table::tstrsplit(age, "\\.")[[1]]
        return(age)
    }
    dt <- data.table::rbindlist(lapply(list.files(path, pattern = type),
                            function(x) fread_agegroup(path, x)))
    return(dt)
}

make_map <- function(dt, type = c("desir", "net")) {
    type <- match.arg(type)
    age <- dt[, unique(agegroup)]
    if(type == "desir") {
        map <- ggplot(set_geom(dt, F)) +
            geom_sf(aes(fill = desirability)) +
            ggtitle(sprintf("%s, estimated desirability", age))
    }
    if(type == "net") {
        map <- ggplot(set_geom(dt, F)) +
            geom_sf(aes(fill = net / pop_all * 100)) +
            ggtitle(sprintf("%s, Net migration", age))
    }
    return(map)
}

make_agemaps <- function(dt) {
    desir_map <- make_map(dt, "desir")
    net_map <- make_map(dt, "net")
    return(desir_map + net_map)
}

dt_geog <- read_output(path, "geog")
dt_dens <- read_output(path, "densfun")
dt_flows <- read_output(path, "flows")
dt_params <- read_output(path, "params")

shp <- setDT(sf::read_sf("~/Diss/inst/extdata/clean/shapes/districts.shp"))
shp[, AGS := as.integer(AGS)]
districts <- fread("./data/districts.csv")

net <- calculate_net(dt_flows, "flows", by = c("year", "agegroup"))
dt_geog[net, net := i.net, on = .(distcode = region, year, agegroup)]
dt_geog[shp, geometry := i.geometry, on = .(distcode = AGS)]
dt_geog[districts, pop_all := i.pop, on = .(distcode, year)]

make_agemaps(dt_geog[agegroup == "below18"]) / make_agemaps(dt_geog[agegroup == "18-25"])
