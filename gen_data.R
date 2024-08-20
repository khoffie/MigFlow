library(data.table)
library(MigStat)
library(sfheaders)
library(sf)

source("functions.R")

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
flows[age_for, frompop := i.german, on = .(fromdist = region, year, agegroup)]
flows[age_for, topop := i.german, on = .(todist = region, year, agegroup)]

shp[, pos := sf::st_point_on_surface(geometry)]
shp[, x := st_coordinates(pos)[, 1]]
shp[, y := st_coordinates(pos)[, 2]]
dt_coords[density, density := i.density, on = .(distcode = region, year)]
dt_coords[age_for[agegroup == "all", ], pop := i.german, on = .(distcode = region, year)]
## dt_coords <- shp[, .(distcode = AGS, year, name = GEN, pop, density, xcoord = x / 1e3, ycoord = y / 1e3,
##                      bl_name, bl_ags)]

### check if districts are the same
all(flows[, unique(fromdist)] == flows[order(todist), unique(todist)] )
all(flows[, unique(fromdist)] == dt_coords[, distcode])

dt_coords[distcode %in% c(5315, 2000, 11000, 14713, 9162, 1001),
          .(name, xcoord = xcoord, ycoord = ycoord)]

regions  <- dt_coords[, unique(distcode)]
distances <- CJ(fromdist = regions, todist = regions)
distances[dt_coords, c("x_from", "y_from") := .(xcoord, ycoord), on = .(fromdist = distcode)]
distances[dt_coords, c("x_to", "y_to") := .(xcoord, ycoord), on = .(todist = distcode)]
distances[, distance := sqrt((x_from - x_to)^2 + (y_from - y_to)^2)]
flows[distances, dist := as.integer(round(i.distance)), on = .(fromdist, todist)]

dt_coords[, .(min = min(xcoord), max = max(xcoord))]
dt_coords[, .(min = min(ycoord), max = max(ycoord))]

dt_coords <- dt_coords[, .(distcode, year, pop, density, name, xcoord, ycoord, bl_ags, bl_name)]
fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")
fwrite(dt_coords, "~/Documents/GermanMigration/data/districts.csv")
