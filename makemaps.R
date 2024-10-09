library(data.table)
library(sf)
library(MigStat)
library(ggplot2)
library(patchwork)
source("functions.R")

path <- "./manuscript_input/bu"


read_agegroup <- function(path, file) {
    get_age <- function(file) {
        age <- data.table::tstrsplit(file, "_")[[2]]
        age <- data.table::tstrsplit(age, "\\.")[[1]]
        return(age)
    }
    dt <- data.table::fread(file.path(path, file))
    dt[, agegroup := get_age(file)]
    return(dt)
}

read_output <- function(path, type) {
    dt <- rbindlist(lapply(list.files(path, pattern = type),
                            function(x) fread_agegroup(path, x)))
    return(dt)
}

dt_geog <- read_output(path, "geog")
dt_dens <- read_output(path, "densfun")
dt_flows <- read_output(path, "flows")
dt_params <- read_output(path, "params")

flows <- data.table::fread(file.path("data", "FlowDataGermans.csv"))
shp <- setDT(sf::read_sf("~/Diss/inst/extdata/clean/shapes/districts.shp"))
districts <- fread("./data/districts.csv")

shp[, AGS := as.integer(AGS)]

net <- calculate_net(flows, "flows", by = c("year", "agegroup"))
net <- net[agegroup == "18-25"]

dt[net, net := i.net, on = .(distcode = region, year, agegroup)]
dt[shp, geometry := i.geometry, on = .(distcode = AGS)]
dt[districts, pop_all := i.pop, on = .(distcode, year)]

desir_map <- ggplot(set_geom(dt, F)) +
    geom_sf(aes(fill = desirability)) +
    ggtitle(sprintf("%s, estimated desirability", age))

net_map <- ggplot(set_geom(dt, F)) +
    geom_sf(aes(fill = net / pop_all * 100)) +
    ggtitle(sprintf("%s, Net migration", age))

desir_map + net_map
