library(data.table)
library(MigStat)

age_for <- fread(file.path(ps$clean, "aux_data", "age17for.csv"))
flows <- ana$flows[year == 2017, .(fromdist = origin, todist = destination,
                                   year, agegroup = age_group, dist = distance_pos,
                                   flows = flow, frompop = pop_o / 1e3, topop = pop_d / 1e3)]
flows[age_for, frompop_ger := i.german / 1e3, on = .(fromdist = region, year, agegroup = age_group)]
flows[age_for, topop_ger := i.german / 1e3, on = .(todist = region, year, agegroup = age_group)]

flows[agegroup == "über65", agegroup := "above65"]
flows[agegroup == "unter18", agegroup := "below18"]
flows <- flows[agegroup != "all"]
fwrite(flows, "~/Documents/GermanMigration/data/FlowDataGermans.csv")
