# MigFlow


## Getting Julia

If you're on Linux, run

```
$ curl -fsSL https://install.julialang.org | sh
```

For Windows and further advice please look
[here](https://github.com/JuliaLang/juliaup). This repository uses
Julia 1.10.4; to get it and make it the default, run


```
$ juliaup add 1.10.4
$ juliaup default 1.10.4
```

Whenever you type `'julia'` into a terminal now, it will start Julia
version 1.10.4.

## Getting all needed libraries
Once you have the correct Julia version, you can clone the repository
and install all needed Julia libraries.

```
$ git clone https://github.com/khoffie/MigFlow.git
$ cd MigFlow
$ julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This will install all libraries needed in the correct version.

## Basic usage

```
using CSV, DataFrames, Turing, StatsBase, Random, Plots, StatsPlots
using ApproxFun, CategoricalArrays, NamedArrays, LaTeXStrings
using ADTypes, ReverseDiff

include("../src/utils.jl")
include("../src/estimation.jl")
include("../src/loadgermdata.jl")
include("../src/diag.jl")
include("../src/diagplots.jl")
include("../src/chebies.jl")
include("../src/gen_mdat.jl")
## available models
include("../src/choiceset.jl")
include("../src/norm.jl")
include("../src/fullmodel.jl")
include("../src/fullmodel2.jl")
include("../src/othermodels.jl")

p = .1 # fraction of rows
data = load_data("30-50", 2017, 1.0, "../data/"; only_positive = true);

## type = "joint" models the joint choice of leaving origin and
## choosing destination, type = "conditional" takes outflux as frompop
## and thus models the conditional choice of destination given origin
## was left.
## ndc "Number of density cheby coefs",
## ngc = "Number of geo cheby coefs"
mdat = gen_mdat(data; type = "joint", distscale = 100.0, ndc = 28, ngc = 1);
out1 = @time estimate(norm, mdat); # MLE

```
