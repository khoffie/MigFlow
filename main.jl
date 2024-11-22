using Revise
includet("src/datafitting.jl")
using RCall
includet("src/postprocess.jl")

settings = Dict(
    :sample_rows => false, # if true 10% sample of rows is used
    :positive_only => true,
    :sampler => externalsampler(SliceSampling.HitAndRun(SliceSteppingOut(2.))),
    ## MH(.1^2*I(78)),
    :sample_size => 10,
    :nchains => 4,
    :thinning => 1,
    :run_optim => false,
    :commit_hash => LibGit2.head("."),
    :fit_us => false,
    :fit_germ => true,
    :distance_type => "pos", ## pos / centroid. pos uses
                            ## sf::st_point_on_surface to calculate
                            ## (xcoord, ycoord). centroid uses
                            ## sf::st_centroid for that. distances
                            ## reflect that
    :topop_type => "agegroup", ## agegroup for age specific population, all for total population
    :year_min => 2017, ## for German data
    :year_max => 2017,
    :agegroups => ["30-50"],
    :outpath => "tempered"
)

function main(settings)
    outpath = makeoutpath(settings[:outpath])
    ## install helpeR only if newer version in repo, never upgrade dependencies
    R"devtools::install_local('./helpeR/', upgrade = 'never', force = FALSE)"
    R"helpeR::gen_data(year_min = $(settings[:year_min]),
                       year_max = $(settings[:year_max]),
                       dist_type = $(settings[:distance_type]),
                       topop_type = $(settings[:topop_type]))"
    mainfit(settings, outpath)
    postprocess(outpath, true)
end

