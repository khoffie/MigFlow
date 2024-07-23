library(data.table)
library(stringr)
library(MigStat)
dt <- fread("~/Documents/GermanMigration/data/opti_vals.csv")
shp <- setDT(sf::read_sf(file.path(ps$clean_shapes, "districts_ext.shp")))
bad_ags <- c(3159, 5978) ## not in julia data
shp <- shp[AGS %notin% bad_ags]
shp[, id := 1:nrow(shp)]


test <- dt[startsWith(parnames, "desir2")]
test[, id := as.numeric(str_extract_all(gsub("desir2", "",parnames), "\\d+"))]

real <- ana$dt_tot[year == 2017 & age_group == "18-25", .(region, real = bal_rel / 10)]
real[, id := 1:399]

test[shp, geometry := i.geometry, on = .(id)]
test[real, real := i.real, on = .(id)]

dt_plt <- melt(test[, .(id, values, real, geometry)], id.vars = c("id", "geometry"))

ggplot(set_geom(dt_plt, F)) +
    geom_sf(aes(fill = value)) +
    facet_wrap(~variable)


dt_plt[, summary(value), keyby = variable]
