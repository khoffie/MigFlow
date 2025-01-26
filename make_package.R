library(devtools)

devtools::create_package("writeup/reporteR")
devtools::document("writeup/reporteR")
devtools::check("writeup/reporteR")
devtools::install("writeup/reporteR", upgrade = FALSE)
