library(data.table)
library(MigStat)
library(sfheaders)
library(sf)

p_clean <- "~/Diss/inst/extdata/clean/"
age_for <- fread(file.path(p_clean, "aux_data", "age17for.csv"))
age_for[age_group == "all", .(sum(all), sum(german))]
shp <- setDT(sf::read_sf(file.path(p_clean, "/shapes/districts_ext.shp")))
shp[, year := 2017] ## kind of hacky but otherwise join fails
density <- data.table::fread(file.path(p_clean, "aux_data", "density.csv"))
density <- density[, .(region, year, density, bl_ags)]
density[is.na(density)]

flows <- ana$flows[year == 2017, .(fromdist = origin, todist = destination,
                                   year, agegroup = age_group, dist = distance_pos,
                                   flows = flow, frompop = pop_o / 1e3, topop = pop_d / 1e3)]
flows <- flows[order(fromdist, todist, year, agegroup)]
flows[age_for, frompop_ger := i.german / 1e3, on = .(fromdist = region, year, agegroup = age_group)]
flows[age_for, topop_ger := i.german / 1e3, on = .(todist = region, year, agegroup = age_group)]
flows[agegroup == "über65", agegroup := "above65"]
flows[agegroup == "unter18", agegroup := "below18"]
flows <- flows[agegroup != "all"]

dt_coords <- shp[, .(distcode = AGS, year, name = GEN, xcoord = x, ycoord = y)]
dt_coords <- dt_coords[distcode %notin% c(3159, 5978)]
dt_coords[density, density := i.density, on = .(distcode = region, year)]

### check if districts are the same
all(flows[, unique(fromdist)] == flows[order(todist), unique(todist)] )
all(flows[, unique(fromdist)] == dt_coords[, distcode])

dt_coords[distcode %in% c(5315, 2000, 11000, 14713, 9162, 1001),
          .(name, xcoord = xcoord / 1000, ycoord = ycoord / 1000)]

fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")
fwrite(dt_coords, "~/Documents/GermanMigration/data/districts.csv")


lapply(flows, function(x) all(is.na(x)))
lapply(dt_coords, function(x) all(is.na(x)))

dt_coords[is.na(density)]

dt_coords
### make sure coords are alright
## need to join name and geometry
## dt_coords[, lbl := substring(name, 1, 2)]
## ggplot(set_geom(dt_coords, F)) +
##     geom_sf() +
##     geom_point(aes(xcoord, ycoord)) +
##     geom_label(aes(xcoord, ycoord, label = lbl))


