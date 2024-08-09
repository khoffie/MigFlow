library(data.table)
library(MigStat)
library(sfheaders)
library(sf)

age_for <- fread(file.path(ps$clean, "aux_data", "age17for.csv"))
shp <- setDT(sf::read_sf(file.path(ps$clean_shapes, "districts_ext.shp")))

flows <- ana$flows[year == 2017, .(fromdist = origin, todist = destination,
                                   year, agegroup = age_group, dist = distance_pos,
                                   flows = flow, frompop = pop_o / 1e3, topop = pop_d / 1e3)]
flows[age_for, frompop_ger := i.german / 1e3, on = .(fromdist = region, year, agegroup = age_group)]
flows[age_for, topop_ger := i.german / 1e3, on = .(todist = region, year, agegroup = age_group)]
## flows[shp, c("x_o", "y_o") := .(i.x, i.y), on = .(fromdist = AGS)]
## flows[shp, c("x_d", "y_d") := .(i.x, i.y), on = .(todist = AGS)]

flows[agegroup == "über65", agegroup := "above65"]
flows[agegroup == "unter18", agegroup := "below18"]
flows <- flows[agegroup != "all"]
fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")
fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")

test <- fread("~/Documents/GermanMigration/data/FlowDataGermans.csv")


dt_coords <- shp[, .(district = AGS, year, xcoord = x, ycoord = y)][district %notin% c(3159, 5978)]
fwrite(dt_coords, "~/Documents/GermanMigration/data/district_coords.csv")
