age(df, age::Vector{AbstractString}) = filter(:agegroup => n -> n ∈ age, df)
age(df, age::AbstractString) = filter(:agegroup => n -> n == age, df)
year(df, y::Vector{Int64}) = filter(:year => n -> n ∈ y, df)
year(df, y::Int64) = filter(:year => n -> n == y, df)
year(df, xmin::Float64, xmax::Float64) = df[df.year .<= xmax .&& df.year .>= xmin, :]
origin(df, o::Vector{Int64}) = filter(:fromdist => n -> n ∈ o, df)
origin(df, o::Int64) = filter(:fromdist => n -> n == o, df)
destination(df, d) = filter(:todist => n -> n ∈ d, df)
code(df, c) = filter(:distcode => n -> n ∈ c, df)
