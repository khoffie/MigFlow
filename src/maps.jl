using Shapefile

shp = Shapefile.Handle("./data/clean/shapes/districts_ext.shp")
f = "germchain_2017_30-50.csv"
path = "./manuscript_input/30kMH"
chain = deserialize(joinpath(path, f))
geo = CSV.read(joinpath(path, "germgeog_2017_30-50.csv"), DataFrame)

xmin = minimum(geo.xcoord)
xmax = maximum(geo.xcoord)
ymin = minimum(geo.ycoord)
ymax = maximum(geo.ycoord)


function plotdesirability(chain, out)
    function minmaxcoords()
        xmin = 303
        xmax = 909
        ymin = 5268
        ymax = 6072
        return (; xmin, xmax, ymin, ymax)
    end
    c = minmaxcoords()
    params = OrderedDict(zip(string.(chain.value.axes[2]), chain.value.data[end, :, 1]))
    coefs = [k for k in keys(params) if contains(k, "desire")]
    desirfun = Fun(ApproxFun.Chebyshev(c.xmin .. c.xmax) * ApproxFun.Chebyshev(c.ymin .. c.ymax),
                   getindex.(Ref(params), coefs) ./ 10)
    des = [desirfun(x, y) for (x, y) in zip(geo.xcoord, geo.ycoord)]
    R"""
    library(ggplot2)
    library(helpeR)
    desirmap <- function(x) {
        shp <- sf::read_sf("./data/clean/shapes/districts_ext.shp")
        shp$desirability <- x
        p <- ggplot(shp) +
            geom_sf(aes(fill = desirability))
        return(p)
    }
    desirmap($des)
    ggsave($(out))
        """
    return des
end

out = "~/Desktop/plt.pdf"

des = plotdesirability(chain, out)
