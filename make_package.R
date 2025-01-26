library(devtools)

devtools::create_package("writeup/reporteR")
devtools::document("reporteR")
devtools::check("reporteR")
devtools::install("reporteR", upgrade = FALSE)
