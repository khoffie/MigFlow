library(helpeR) ## devtools::install_local("helpeR"), installs other packages as well

## probably a good idea to add directories automatically
## raw/german_flows
## clean/german_flows


rawdata_to_cleandata <- function(raw, clean) {
  dir.create(file.path(clean))
  dir.create(file.path(clean, "german_flows"))
  helpeR::preprocess_germanflows(file.path(raw, "german_flows"),
                                 file.path(clean, "german_flows"))
  ## This may create an NA.csv which makes german_flows() fail.
  helpeR::german_flows(file.path(clean, "german_flows"),
                       file.path(clean, "flows_districts_2000_2017_ger.csv"))

  helpeR::german_popdata(raw, clean)

  ## German auxilliary data
  german_administrative_changes(file.path(raw, "ref-kreise-umrech-2019-1990-2018.xlsx"),
                                file.path(clean, "correct.csv"))
  file.copy(file.path(raw, "density.csv"), file.path(clean, "density.csv"),
            overwrite = TRUE)
  file.copy(file.path(raw, "shapes"), clean,
            recursive = TRUE, overwrite = TRUE)
}

## data preprocessing
rawdata_to_cleandata("./data/raw", "./data/clean")

## generates data used for modeling: FlowDataGermans.csv and districts.csv
helpeR::gen_data("./data", year_min =  2000, year_max = 2017, dist_type = "pos", topop_type = "all")
