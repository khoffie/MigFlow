# using Pkg
# Pkg.activate(".")
using Revise
includet("src/datafitting.jl")
using RCall
includet("src/postprocess.jl")

settings = Dict(
    :model_type => "base", # base: Fundamentals model, full:
                          # fundamentals + density + desirability
    :sample_rows => false, # if true 10% sample of rows is used
    :positive_only => true,
    :sampler => "MH", # "MH" or "Slice"
    # MH(.1^2*I(6)),
    ## externalsampler(SliceSampling.HitAndRun(SliceSteppingOut(2.))),
    :sample_size => 10,
    :nchains => 4,
    :thinning => 1,
    ## :run_optim => false,
    :commit_hash => LibGit2.head("."),
    :fit_us => false,
    :fit_germ => true,
    :distance_type => "pos", ## pos / centroid. pos uses
                            ## sf::st_point_on_surface to calculate
                            ## (xcoord, ycoord). centroid uses
                            ## sf::st_centroid for that. distances
                            ## reflect that
    :topop_type => "all", ## agegroup for age specific population, all for total population
    :year_min => 2017, ## for German data
    :year_max => 2017,
    :agegroups => nothing, # nothing: all age groups
    :outpath => "fundamentals",
    :rm_dir => true,
    :temp_samples => 10,
    :min_temp => 8500,
    :temp_decay => 0.2
)

function main(settings)
    outpath = makeoutpath(settings[:outpath])
    ## install helpeR only if newer version in repo, never upgrade dependencies
    ## R"devtools::install_local('./helpeR/', upgrade = 'never', force = FALSE)"
    # R"helpeR::gen_data(year_min = $(settings[:year_min]),
    #                    year_max = $(settings[:year_max]),
    #                    dist_type = $(settings[:distance_type]),
    #                    topop_type = $(settings[:topop_type]))"

    run(`R -q -e 'helpeR::gen_data(
    year_min = 2017,
    year_max = 2017,
    dist_type = "pos",
    topop_type = "all")'`)

    mainfit(settings, outpath)
    postprocess(path = outpath, render_doc = false)
end

main(settings)

# R"helpeR::gen_data(year_min = 2017,
#                        year_max = 2017,
#                        dist_type = 'pos',
#                        topop_type = 'all')"


# R -q -e 'helpeR::gen_data(year_min = 2017, year_max = 2017, dist_type = "pos", topop_type = "all")'
