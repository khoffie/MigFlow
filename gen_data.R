library(data.table)
library(MigStat)

age_for <- fread(file.path(ps$clean, "aux_data", "age17for.csv"))
flows <- ana$flows[year == 2017, .(fromdist = origin, todist = destination,
                                   year, agegroup = age_group, dist = distance_pos,
                                   flows = flow, frompop = pop_o, topop = pop_d)]
flows[age_for, frompop_ger := i.german, on = .(fromdist = region, year, agegroup = age_group)]
flows[age_for, topop_ger := i.german, on = .(todist = region, year, agegroup = age_group)]

fwrite(flows, "~/Documents/GermanMigration/FlowDataGermans.csv")
