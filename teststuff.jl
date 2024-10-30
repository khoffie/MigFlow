using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots
germ = loadallGermData(; sample = false)
us = loadallUSdata()
