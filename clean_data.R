library(data.table)
library(helpeR)
library(sf)
library(readxl)

helpeR::german_popdata("./data/raw", "./data/clean")

raw <- file.path("./data/raw", "german_flows")
clean <- file.path("./data/clean", "german_flows")

## shp <- data.table(sf::read_sf("./data/clean/shapes/districts_ext.shp"))

flows <- rbindlist(lapply(file.path(clean, list.files(clean)), fread))
flows[year != 2003] ## there are lots of problems

## I think every unique combination of keys should be in the data only
## once. This seems to be wrong
keys <- c("origin", "destination", "age_group", "year")
flows[, .N, keyby = keys][N > 1]
flows[, sum(flow), keyby = keys]
flows[, .N, keyby = c(keys, "flow")][N > 1][, unique(year)]
flows <- flows[, .(flow = sum(flow)), keyby = keys]
## Just sum over them. This assumes all rows, even the dublicated and
## even the ones where even the flow is dublicated, are valid. What
## difference does this make compared to taking only one
##flows[, .SD[1], keyby = keys][, sum(flow)]

flows <- flows[year != 2003, ] ## problems with age_groups in BadenWürrtemberg2003
## there are some ags to fix. From 2000Hamburg I discovered that there
## are different AGS for every quarter of Hamburg, I remap those
ags_hamburg <- 2101:2718
flows[origin %in% ags_hamburg, "origin" := 2000]
flows[destination %in% ags_hamburg, "destination" := 2000]
ags_berlin <- c(11001:11012, 11100:11223)
flows[origin %in% ags_berlin, "origin" := 11000]
flows[destination %in% ags_berlin, "destination" := 11000]
years <- 2000 : 2017
lapply(years, function(y) setdiff(flows[year == y, unique(origin)], shp[year == 2018, AGS]))
lapply(years, function(y) setdiff(flows[year == y, unique(destination)], shp[year == 2018, AGS]))
#### now many ags that were different before are the same now so I sum
#### again by keys
flows <- flows[, .(flow = sum(flow)), keyby = keys]
### age groups are often named differently


## remove weird age groups
flows <- flows[age_group != "AG", ] ## dont know whats up here
## only in 2006, seems to be more or less the sum of the other age groups
flows <- flows[age_group != "insgesamt", ] 
flows[age_group == "18-25", "age_group" := "18 - 25"]
flows[age_group == "25-30", "age_group" := "25 - 30"]
flows[age_group == "30-50", "age_group" := "30 - 50"]
flows[age_group == "50-65", "age_group" := "50 - 65"]
flows[age_group == "18 bis unter 25", "age_group" := "18 - 25"]
flows[age_group == "25 bis unter 30", "age_group" := "25 - 30"]
flows[age_group == "30 bis unter 50", "age_group" := "30 - 50"]
flows[age_group == "50 bis unter 65", "age_group" := "50 - 65"]
flows[age_group == "65 und mehr", "age_group" := "65 und älter"]
flows[age_group == "65 und älter", "age_group" := "über 65"]
flows[, uniqueN(age_group)] == 6

### It is useful to have an age group consisting of all age groups
dt_all_ages <- flows[, .(flow = sum(flow)), keyby = c("origin", "destination", "year")]
dt_all_ages[, "age_group" := "all"]
keys <- c("origin", "destination", "year", "age_group")
setkeyv(flows, keys)
setkeyv(dt_all_ages, keys)
keys_all <- rbind(flows[, ..keys], dt_all_ages[, ..keys])
flows <- flows[keys_all]
setkeyv(flows, keys)
## should be doable easier
flows[dt_all_ages, "flow":= i.flow]

if(flows[age_group != "all", sum(flow)] != flows[age_group == "all", sum(flow)]) {
    stop("Age group all not sum of others!")
}

flows[, "age_group" := gsub(" ", "", age_group)] ## do in preparation
## I confused origin and destination when preparing the data. Actually
## not I confused them but they are confused in the original data
setnames(flows, c("origin", "destination"), c("destination", "origin"))
setcolorder(flows, c("origin", "destination"))

## removes flows within Berlin and Hamburg since in the data they have
## separate AGS for their quarters. Also rows are removed where origin
## == destination but 0 flows. Not sure why they are in the data in
## the first place. Maybe because of age_group stuff?
flows <- flows[origin != destination]
