library(data.table)
library(MigStat)
library(sfheaders)
library(sf)

p_clean <- "~/Diss/inst/extdata/clean/"
flows <- fread(file.path(p_clean, "flows_districts/districts_2000_2017_ger.csv"))
age_for <- fread(file.path(p_clean, "aux_data", "age17for.csv"))
age_for[age_group == "all", .(sum(all), sum(german))]
shp <- setDT(sf::read_sf(file.path(p_clean, "/shapes/districts_ext.shp")))
shp[, year := 2017] ## kind of hacky but otherwise join fails
density <- data.table::fread(file.path(p_clean, "aux_data", "density.csv"))
density <- density[, .(region, year, density, bl_ags)]

flows <- flows[year == 2017, .(fromdist = origin, todist = destination,
                                   year, agegroup = age_group, flows = flow)]

flows <- flows[order(fromdist, todist, year, agegroup)]
flows[age_for, frompop_ger := i.german / 1e3, on = .(fromdist = region, year, agegroup = age_group)]
flows[age_for, topop_ger := i.german / 1e3, on = .(todist = region, year, agegroup = age_group)]
flows[agegroup == "über65", agegroup := "above65"]
flows[agegroup == "unter18", agegroup := "below18"]
flows <- flows[agegroup != "all"]

shp[, pos := sf::st_point_on_surface(geometry)]
shp[, x := st_coordinates(pos)[, 1]]
shp[, y := st_coordinates(pos)[, 2]]
dt_coords <- shp[, .(distcode = AGS, year, name = GEN, xcoord = x / 1e3, ycoord = y / 1e3)]
dt_coords[density, density := i.density, on = .(distcode = region, year)]

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
flows[distances, distance := as.integer(round(i.distance)), on = .(fromdist, todist)]

dt_coords
dt_coords[, .(min = min(xcoord), max = max(xcoord))]
dt_coords[, .(min = min(ycoord), max = max(ycoord))]

fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")
fwrite(dt_coords, "~/Documents/GermanMigration/data/districts.csv")

