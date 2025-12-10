using CSV, DataFrames

include("../src/utils.jl")
df = CSV.read("../data/FlowDataGermans.csv", DataFrame)
districts = CSV.read("../data/districts.csv", DataFrame)

function archive_data(df, di, path, name)
    p = joinpath(path, "clean/shapes")
    if ! isdir(p); mkpath(p); end
    CSV.write(joinpath(path, "FlowDataGermans.csv"), df)
    CSV.write(joinpath(path, "districts.csv"), di)
    for name in ["districts_ext", "states"], ext in ["dbf", "prj", "shp", "shx"]
        src = joinpath("../data/clean/shapes", "$name.$ext")
        dst = joinpath(p, "$name.$ext")
        cp(src, dst, force = true)
    end
    run(`tar -czvf $(name) -C $(path) .`)
    rm("../data/archives", recursive = true)
end

archive_data(age(year(df, 2017), "18-25"), year(districts, 2017),
             "../data/archives/minimal", "../data/data_minimal.tar.gz")

archive_data(df, districts, "../data/archives/full", "../data/data_full.tar.gz")

run(`tar -czvf ../data/data_raw.tar.gz -C ../data/raw .`)
