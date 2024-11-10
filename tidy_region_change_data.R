library("data.table")
library("readxl")

#################### Data and Functions ##############################
## p <- "~/Diss/exec/region_changes/"
## p <- "~/Diss/exec/one_run//data/raw/corrections"
p <- paths$raw_corrections
changes_f <- "ref-kreise-umrech-2019-1990-2018.xlsx"
changes_f <- "districts_19.csv"
p_out <- "~/Diss/exec/one_run//data/clean/corrections"
p_out <- paths$clean_corrections
p_shps <- "~/Diss/exec/one_run/data/clean/shapes/"
p_shps <- paths$clean_shapes
shp <- MigStat:::read_clean_shps(p_shps, "complete")$districts
excel_sheet_reader <- function(filename) {
  sheets <- excel_sheets(filename)
  x <- lapply(sheets, function(X) read_excel(filename, sheet = X))
  names(x) <- sheets
  return(x)
}

################### Making one table #################################
changes <- excel_sheet_reader(file.path(p, changes_f))
## colnames(changes[["2000-2019"]])[2] <- "name_old"
## changes[["2000-2019"]][name_old == "KS Berlin"]
## okay so read_xl seems to read in the ags wrongly
lapply(changes, setDT)
for (i in seq_along(changes)) {
    changes[[i]][, "year" := names(changes)[i]]
}
## before I can comnine the lists to one table I need to make sure
## they have the same column names
for (i in seq_along(changes)) {
    dt <- changes[[i]]
    setnames(dt, colnames(dt)[1], "ags_old")
    setnames(dt, colnames(dt)[2], "name_old")
    old <- "flächen-\r\nproportionaler\r\nUmsteige-schlüssel"
    setnames(dt, old, "conv_f")
    old <- "bevölkerungs- \r\nproportionaler \r\nUmsteige- \r\nschlüssel"
    setnames(dt, old, "conv_p")
    old <- "Kreise\r\n 31.12.2019"
    setnames(dt, old, "ags_new")
    old <- "Kreisname 2019"
    setnames(dt, old, "name_new")
}
keep <- c("ags_old", "name_old", "conv_f", "conv_p",
          "ags_new", "name_new", "year")
for (i in seq_along(changes)) {
    changes[[i]] <- changes[[i]][, ..keep]
}
changes <- rbindlist(changes)
### year can be better handled if -2019 is omitted
changes[, "year" := tstrsplit(year, "-", keep = 1)]
### some manual corrections are necessary, might be because read_excel
### reads some ags wrongly
changes[, "ags_old" := as.integer(gsub("000$", "", as.character(ags_old)))]
changes[, "ags_new" := as.integer(gsub("000$", "", as.character(ags_new)))]
changes[ags_old == 2000000, "ags_old" := 2000]
changes[ags_new == 2000000, "ags_new" := 2000]
changes[ags_old == 11000000, "ags_old" := 11000]
changes[ags_new == 11000000, "ags_new" := 11000]

### check if all ags can be found in the shapefiles to make sure all
### data us the same format
## lapply(2000:2018, function(y) setdiff(shp[year == y, unique(AGS)],
##                                       changes[year == y, unique(ags_old)]))
not_found <- lapply(2000:2018, function(y)
    length(setdiff(changes[year == y, unique(ags_old)],
                   shp[year == y, unique(AGS)])))

## better to only save when really all ags are found
if (Reduce(sum, not_found) != 0) {
    stop("Some AGS not found in shapefiles")
}
if (Reduce(sum, not_found) == 0) {
    data.table::fwrite(changes, file.path(p_out, changes_f))
}

