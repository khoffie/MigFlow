using StatsBase
germ = loadallGermData(; sample = false)
us = loadallUSdata()


sample = true    
test = loadallGermData(0; sample = sample)
test = loadallUSdata(0; sample = sample)
nzeros = 0
names(geog)

@which maximum_a_posteriori

?Turing.Optimisation.maximum_a_posteriori
