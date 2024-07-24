## A file for Julia functions that take Julia inputs and call R scripts to plot or analyze them
## first we load libraries, then we make some functions

R"""
library(ggplot2) 

"""

function plotdesirability(desir)

    @rput desir

# do plotting here in R
R"""

plot(desir$x,netmigration$y)

"""
end
