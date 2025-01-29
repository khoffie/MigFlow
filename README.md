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
