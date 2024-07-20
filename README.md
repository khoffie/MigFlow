# GermanMigration
Study of German Migration patterns

## Getting Julia

To get started download Julia from the julialang site and extract it in your home directory, and then link your home dir's bin/julia to the downloaded version. Assuming you have ~/bin/ in your path, you can then run julia easily.

```
$ wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.4-linux-x86_64.tar.gz
$ tar zxf julia-1.10.4-linux-x86_64.tar.gz
$ ln -s ~/julia-1.10.4/bin/julia bin/julia

```

## Cloning the Repo and linking to your data file

replace <path_to_data_file> with a path to a data file set up like what's described in fitturing.jl in the 
comments 

```
git clone https://github.com/dlakelan/GermanMigration.git
cd GermanMigration
ln -s <path_to_data_file> data/FlowData.csv
```

## Setting up the Julia environment

To set up the environment after cloning the repo from inside the GermanMigration directory at the shell prompt do:

`$ julia ./setupenv.jl`

This will download all the packages

## Now to explore the data:

```
julia

...
julia> include("fitturing.jl")
```

From here you can explore the two fit models or debug whatever bugs are in my script :-)

