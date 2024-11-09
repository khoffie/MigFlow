library(data.table)
library(helpeR)
library(sf)
library(readxl)


helpeR::preprocess_germanflows("./data/raw/german_flows", "./data/clean/german_flows")
helpeR::german_flows("./data/clean/german_flows", "./data/clean")

helpeR::german_popdata("./data/raw", "./data/clean")
