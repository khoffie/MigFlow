library(data.table)
library(helpeR)
library(sfheaders)
library(sf)

read_age <- function() {
    dt <- fread(file.path(p_clean, "aux_data", "age17for.csv"))
    data.table::setnames(dt, "age_group", "agegroup")
    helpeR::rec_ages(dt)
    return(dt)
}

p_clean <- "~/Diss/inst/extdata/clean/"
flows <- fread(file.path(p_clean, "flows_districts/districts_2000_2017_ger.csv"))
age_for <- read_age()
shp <- setDT(sf::read_sf(file.path(p_clean, "/shapes/districts_ext.shp")))
density <- data.table::fread(file.path(p_clean, "aux_data", "density.csv"))[
                         , .(region, year, density, bl_ags)]

clean_flows <- function(flows, age_dt) {
    flows <- flows[year == 2017 & age_group != "all", 
               .(fromdist = origin, todist = destination, year,
                 agegroup = age_group, flows = flow)]
    rec_ages(flows)
    flows <- helpeR::add_missing_flows(flows, flows[, unique(fromdist)],
                              flows[, unique(agegroup)], flows[, unique(year)])
    ### omitting since all intra-district flows are 0 and this would be
    ### hard to explain for the model
    flows <- flows[fromdist != todist]
    flows <- flows[order(fromdist, todist, year, agegroup)]
    
    flows[age_dt, frompop := i.german, on = .(fromdist = region, year, agegroup)]
    flows[age_dt, topop := i.german, on = .(todist = region, year, agegroup)]
    return(flows)
}

gen_coords <- function(dt, type = c("centroid", "pos")) {
    t <- match.arg(type)
    if(t == "centroid") {
        dt[, pos := sf::st_centroid(geometry)]
    }
    if(t == "pos") {
        dt[, pos := sf::st_point_on_surface(geometry)]
    }
    dt[, xcoord := st_coordinates(pos)[, 1] / 1e3]
    dt[, ycoord := st_coordinates(pos)[, 2] / 1e3]
    return(dt)
}

gen_coords_dt <- function(shp, age_dt, density_dt, type) {
    shp[, year := 2017] ## actual year is 2018, hacky but otherwise join fails
    gen_coords(shp, type = type)
    dt <- shp[, .(distcode = AGS, year, name = GEN, 
                         xcoord, ycoord,
                         bl_name, bl_ags)]

    dt[density, density := i.density, on = .(distcode = region, year)]
    dt[age_for[agegroup == "all", ], pop := i.german, on = .(distcode = region, year)]
    dt <- dt[, .(distcode, year, pop, density, name, xcoord, ycoord, bl_ags, bl_name)]
    return(dt)
}

check_tables <- function(flows, coords) {
    stopifnot(all(flows[, unique(fromdist)] == flows[order(todist), unique(todist)] ))
    stopifnot(all(flows[, unique(fromdist)] == coords[, distcode]))
    coords[distcode %in% c(5315, 2000, 11000, 14713, 9162, 1001),
              .(name, xcoord = xcoord, ycoord = ycoord)]
    coords[, .(min = min(xcoord), max = max(xcoord))]
    coords[, .(min = min(ycoord), max = max(ycoord))]
}

calculate_distances <- function(flows, coords) {
    regions  <- coords[, unique(distcode)]
    distances <- CJ(fromdist = regions, todist = regions)
    distances[coords, c("x_from", "y_from") := .(xcoord, ycoord), on = .(fromdist = distcode)]
    distances[coords, c("x_to", "y_to") := .(xcoord, ycoord), on = .(todist = distcode)]
    distances[, distance := sqrt((x_from - x_to)^2 + (y_from - y_to)^2)]
    flows[distances, dist := as.integer(round(i.distance)), on = .(fromdist, todist)]
    return(NULL)
}

pop_weighted_distance <- function() {
    gen_munis <- function() {
        munis <- setDT(sf::read_sf("~/Diss/inst/extdata/clean/shapes/munis_all.shp"))[year == 2017]
        ink <- fread("~/Diss/inst/extdata/clean/inkar/inkar_2021.csv")
        munis_pop <- modeleR::create_design_mat(ink, 425, "Gemeinden", 2017)
        munis_pop[X425 == 0.1, X425 := 0] ## because create_design_mat sets them to 0.1
        munis[, distcode := as.integer(substr(AGS, 1, 5))]
        munis[, AGS := as.integer(AGS)]
        munis[munis_pop, pop := i.X425, on = .(AGS)]
        gen_coords(munis)
        munis[, pop_dist := sum(pop, na.rm = TRUE), keyby = .(distcode)]
        munis <- munis[, .(municode = AGS, distcode, pop, pop_dist, xcoord, ycoord)]
        munis <- munis[!is.na(pop)]
        return(munis)
    }

    gen_dist_muni_dt <- function(munis) {
        distances <- data.table::CJ(fromdist = districts[, distcode], todist = districts[, distcode])
        distances <- distances[fromdist != todist]
        cols <- setdiff(colnames(munis), "distcode")
        distances <- distances[munis, on = .(fromdist = distcode), allow.cartesian = TRUE]
        setnames(distances, cols, paste0("from_", cols))
        distances <- distances[munis, on = .(todist = distcode), allow.cartesian = TRUE]
        setnames(distances, cols, paste0("to_", cols))
        return(distances)
    }

    calc_weighted_distance <- function(dt) {
         dt[, distance := sqrt((from_xcoord - to_xcoord)^2 + (from_ycoord - to_ycoord)^2)]
         dt[, fromweights := from_pop / from_pop_dist]
         dt[, toweights := to_pop / to_pop_dist]
         dt <- dt[, .(distance_w = sum(distance * fromweights * toweights)), keyby = .(fromdist, todist)]
         return(dt)
    }
    munis <- gen_munis()
    distances <- gen_dist_muni_dt(munis)
    distances <- calc_weighted_distance(distances)
    return(distances)
}

flows <- clean_flows(flows, age_for)
districts <- gen_coords_dt(shp, age_for, density, type = "pos")
check_tables(flows, districts)
calculate_distances(flows, districts)
## distances <- pop_weighted_distance() ## would be good to calc all distances here
## flows[distances, dist_w := i.distance_w, on = .(fromdist, todist)]

fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")
fwrite(districts, "~/Documents/GermanMigration/data/districts.csv")
