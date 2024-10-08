library(data.table)
library(MigStat)
library(sfheaders)
library(sf)

source("functions.R")

p_clean <- "~/Diss/inst/extdata/clean/"
flows <- fread(file.path(p_clean, "flows_districts/districts_2000_2017_ger.csv"))
age_for <- fread(file.path(p_clean, "aux_data", "age17for.csv"))
setnames(age_for, "age_group", "agegroup")
rec_ages(age_for)
shp <- setDT(sf::read_sf(file.path(p_clean, "/shapes/districts_ext.shp")))

density <- data.table::fread(file.path(p_clean, "aux_data", "density.csv"))[
                         , .(region, year, density, bl_ags)]

clean_flows <- function(flows, age_dt) {
    flows <- flows[year == 2017 & age_group != "all", 
               .(fromdist = origin, todist = destination, year,
                 agegroup = age_group, flows = flow)]
    rec_ages(flows)
    flows <- add_mising_flows(flows, flows[, unique(fromdist)],
                              flows[, unique(agegroup)], flows[, unique(year)])
### omitting since all intra-district flows are 0 and this would be
### hard to explain for the model
    flows <- flows[fromdist != todist]
    flows <- flows[order(fromdist, todist, year, agegroup)]
    
    flows[age_dt, frompop := i.german, on = .(fromdist = region, year, agegroup)]
    flows[age_dt, topop := i.german, on = .(todist = region, year, agegroup)]
    return(flows)
}

gen_coords_dt <- function(shp, age_dt, density_dt) {
    shp[, year := 2017] ## actual year is 2018, hacky but otherwise join fails
    shp[, pos := sf::st_point_on_surface(geometry)]
    shp[, x := st_coordinates(pos)[, 1]]
    shp[, y := st_coordinates(pos)[, 2]]
    dt <- shp[, .(distcode = AGS, year, name = GEN, 
                         xcoord = x / 1e3, ycoord = y / 1e3,
                         bl_name, bl_ags)]

    dt[density, density := i.density, on = .(distcode = region, year)]
    dt[age_for[agegroup == "all", ], pop := i.german, on = .(distcode = region, year)]
    dt_coords <- dt_coords[, .(distcode, year, pop, density, name, xcoord, ycoord, bl_ags, bl_name)]
    return(dt)
}

check_tables <- function(flows, coords_dt) {
    stopifnot(all(flows[, unique(fromdist)] == flows[order(todist), unique(todist)] ))
    stopifnot(all(flows[, unique(fromdist)] == dt_coords[, distcode]))
    dt_coords[distcode %in% c(5315, 2000, 11000, 14713, 9162, 1001),
              .(name, xcoord = xcoord, ycoord = ycoord)]
    dt_coords[, .(min = min(xcoord), max = max(xcoord))]
    dt_coords[, .(min = min(ycoord), max = max(ycoord))]
}

calculate_distances <- function(flows, coords_dt) {
    regions  <- dt_coords[, unique(distcode)]
    distances <- CJ(fromdist = regions, todist = regions)
    distances[dt_coords, c("x_from", "y_from") := .(xcoord, ycoord), on = .(fromdist = distcode)]
    distances[dt_coords, c("x_to", "y_to") := .(xcoord, ycoord), on = .(todist = distcode)]
    distances[, distance := sqrt((x_from - x_to)^2 + (y_from - y_to)^2)]
    flows[distances, dist := as.integer(round(i.distance)), on = .(fromdist, todist)]
}

flows <- clean_flows(flows, age_for)
dt_coords <- gen_coords_dt(shp, age_for, density)
check_tables(flows, dt_coords)
calculate_distances(flows, coords_dt)

fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")
fwrite(dt_coords, "~/Documents/GermanMigration/data/districts.csv")
