library(helpeR) ## devtools::install_local("helpeR"), installs other packages as well
library(data.table)
library(sf)
library(readxl)

rawdata_to_cleandata <- function(raw, clean) {
    helpeR::preprocess_germanflows(file.path(raw, "german_flows"),
                                   file.path(clean, "german_flows"))
    helpeR::german_flows(file.path(clean, "german_flows"),
                         file.path(clean, "flows_districts_2000_2017_ger.csv"))

    helpeR::german_popdata(raw, clean)

    ## German auxilliary data
    file.copy(file.path(raw, "correct.csv"), file.path(clean, "correct.csv"),
              overwrite = TRUE)
    file.copy(file.path(raw, "density.csv"), file.path(clean, "density.csv"),
              overwrite = TRUE)
    file.copy(file.path(raw, "shapes"), file.path(clean, "shapes"),
              recursive = TRUE, overwrite = TRUE)
}
