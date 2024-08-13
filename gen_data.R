library(data.table)
library(MigStat)
library(sfheaders)
library(sf)

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

p_clean <- "~/Diss/inst/extdata/clean/"
flows <- fread(file.path(p_clean, "flows_districts/districts_2000_2017_ger.csv"))
age_for <- fread(file.path(p_clean, "aux_data", "age17for.csv"))
shp <- setDT(sf::read_sf(file.path(p_clean, "/shapes/districts_ext.shp")))
density <- data.table::fread(file.path(p_clean, "aux_data", "density.csv"))

flows <- flows[year == 2017 & age_group != "all", 
               .(fromdist = origin, todist = destination, year,
                 agegroup = age_group, flows = flow)]

rec_ages(flows)
setnames(age_for, "age_group", "agegroup")
rec_ages(age_for)
shp[, year := 2017] ## actual year is 2018, hacky but otherwise join fails
density <- density[, .(region, year, density, bl_ags)]

flows <- add_mising_flows(flows, flows[, unique(fromdist)],
                          flows[, unique(agegroup)], flows[, unique(year)])
### omitting since all intra-district flows are 0 and this would be
### hard to explain for the model
flows <- flows[fromdist != todist]


flows <- flows[order(fromdist, todist, year, agegroup)]
flows[age_for, frompop := i.german / 1e3, on = .(fromdist = region, year, agegroup)]
flows[age_for, topop := i.german / 1e3, on = .(todist = region, year, agegroup)]

shp[, pos := sf::st_point_on_surface(geometry)]
shp[, x := st_coordinates(pos)[, 1]]
shp[, y := st_coordinates(pos)[, 2]]
dt_coords <- shp[, .(distcode = AGS, year, name = GEN, xcoord = x / 1e3, ycoord = y / 1e3)]
dt_coords[density, density := i.density, on = .(distcode = region, year)]
dt_coords[age_for[agegroup == "all", ], pop := i.german / 1e3, on = .(distcode = region, year)]

### check if districts are the same
all(flows[, unique(fromdist)] == flows[order(todist), unique(todist)] )
all(flows[, unique(fromdist)] == dt_coords[, distcode])

dt_coords[distcode %in% c(5315, 2000, 11000, 14713, 9162, 1001),
          .(name, xcoord = xcoord, ycoord = ycoord)]
### make sure coords are alright
## need to join name and geometry
## dt_coords[, lbl := substring(name, 1, 2)]
## ggplot(set_geom(dt_coords, F)) +
##     geom_sf() +
##     geom_point(aes(xcoord, ycoord)) +
##     geom_label(aes(xcoord, ycoord, label = lbl))

regions  <- dt_coords[, unique(distcode)]
distances <- CJ(fromdist = regions, todist = regions)
distances[dt_coords, c("x_from", "y_from") := .(xcoord, ycoord), on = .(fromdist = distcode)]
distances[dt_coords, c("x_to", "y_to") := .(xcoord, ycoord), on = .(todist = distcode)]
distances[, distance := sqrt((x_from - x_to)^2 + (y_from - y_to)^2)]
flows[distances, dist := as.integer(round(i.distance)), on = .(fromdist, todist)]

dt_coords[, .(min = min(xcoord), max = max(xcoord))]
dt_coords[, .(min = min(ycoord), max = max(ycoord))]

fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")
fwrite(dt_coords, "~/Documents/GermanMigration/data/districts.csv")

library(ggplot2)

ggplot(flows[year == 2017], aes(log(flows))) +
    geom_density() +
    facet_wrap(vars(agegroup)) +
    theme_minimal()

ggplot(flows[year == 2017], aes(log(flows))) +
    geom_histogram() +
    facet_wrap(vars(agegroup)) +
    theme_minimal()

ggplot(flows[year == 2017], aes((flows))) +
    geom_histogram() +
    xlim(c(0, 50)) +
    facet_wrap(vars(agegroup)) +
    theme_minimal()

log(1)
