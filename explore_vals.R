library(data.table)
library(ggplot2)
library(stringr)
library(MigStat)
dt <- fread("~/Documents/GermanMigration/data/opti_vals.csv")
shp <- setDT(sf::read_sf(file.path(ps$clean_shapes, "districts_ext.shp")))
bad_ags <- c(3159, 5978) ## not in julia data
shp <- shp[AGS %notin% bad_ags]
shp[, id := 1:.N]

dt[, param := gsub("desir", "", param)]
dt[, tstrsplit(param, "\\[")]



ggplot(test, aes(estim)) +
    geom_density() +
    facet_wrap(~agegroup)



test[, id := as.numeric(str_extract_all(gsub("desir2", "",param), "\\d+"))]
test <- test[id %notin% very_low]

real <- ana$dt_tot[year == 2017 & age_group == "18-25", .(region, real = bal_rel / 10)]
real[, id := 1:399]


test[shp, geometry := i.geometry, on = .(id)]
test[real, real := i.real, on = .(id)]


test[order(estim)]
shp[id == 257]
very_low <- c(257, 46, 310)
plot(test[, .(estim, real)])

dt_plt <- melt(test[, .(id, estim, real, geometry)], id.vars = c("id", "geometry"))
dt_plt[, value := scale(value), keyby = .(variable)]

plt_estim <- ggplot(set_geom(dt_plt[id %notin% very_low & variable == "estim"], F)) +
    geom_sf(aes(fill = value)) +
    facet_wrap(~variable)

plt_real <- ggplot(set_geom(dt_plt[variable == "real"], F)) +
    geom_sf(aes(fill = value)) +
    facet_wrap(~variable)

plt_real + plt_estim


dt_plt[, .(mean = mean(value), sd = sd(value)), keyby = variable]

    a ~ filldist(Gamma(5.0, 10.0/4.0),Nages)
    b ~ filldist(Gamma(3.0, 1.0/2.0),Nages)
    c ~ filldist(Gamma(5.0, 1.0/4.0),Nages)
    d0 ~ filldist(Gamma(5.0, 0.10/4.0),Nages)
    neterr ~ Gamma(3.0, 0.1/2.0)

    desir ~ filldist(Gamma(400.0, 1.0/399.0), Ndist, Nages)


x <- seq(-100, 100, by = .01)
dt <- data.table(x = x)
dt[, a := dgamma(x, shape = 5, scale = 2.5)]
dt[, b := dgamma(x, shape = 3, scale = .5)]
dt[, c := dgamma(x, shape = 5, scale = .25)]
dt[, d0 := dgamma(x, shape = 5, scale = .025)]
dt[, neterr := dgamma(x, shape = 3, scale = .5)]
dt[, desir := dgamma(x, shape = 400, scale = 1/399)]
dtm <- melt(dt, id.vars = "x")
dtm[value > 0.0001]


ggplot(dtm[value > 0.0001], aes(x, value)) +
    geom_line() +
    facet_wrap(~variable, scale = "free") +
    ggtitle("Prior distributions") +
    theme_minimal()
