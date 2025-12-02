getageyear(r) = r.mdl.meta.age, r.mdl.meta.year
getmodel(r) = r.mdl.meta.model

function addlines!(ax, df, group, xcol, ycol, lw = 1)
    groups = unique(df[!, group])
    for g in groups
        foo = df[df[!, group] .== g, :]
        kwargs = Dict{Symbol, Any}()
        if "col" in names(foo)
            kwargs[:color] = foo.col
        end
        if "lw" in names(foo)
            kwargs[:linewidth] = foo.lw
        end
        lines!(ax, foo[!, xcol], foo[!, ycol], label = g, linewidth = lw; kwargs...)
    end
end

function addlines!(ax, df, group, xcol, ycol, lw = 1)
    return addgroups!(ax, lines!, df, group, xcol, ycol, lw)
end

function addgroups!(ax, f, df, group, xcol, ycol, lw = 1; kwargs...)
    groups = unique(df[!, group])
    for g in groups
        foo = df[df[!, group] .== g, :]
        inner_kwargs = Dict{Symbol, Any}()
        if "col" in names(foo)
            inner_kwargs[:color] = foo.col
        end
        if "lw" in names(foo)
            inner_kwargs[:linewidth] = foo.lw
        end
        kwargs = merge(inner_kwargs, Dict(kwargs...))
        if isnothing(ycol)
            f(ax, foo[!, xcol], label = g; kwargs...)
        else
            f(ax, foo[!, xcol], foo[!, ycol], label = g, linewidth = lw; kwargs...)
        end
    end
end

function grid_position(k::Int, ncols::Int=2)
    i = div(k - 1, ncols) + 1   # row index
    j = mod(k - 1, ncols) + 1   # column index
    return (i, j)
end

function prettytitle!(f, title, subtitle = nothing, row = 0, size = 12)
    pt = 4 / 3
    size = size * pt
    titlelayout = GridLayout(f[row, 1:end], tellwidth = false)
    Label(titlelayout[1, 1], title, fontsize = size,
          font = "TeX Gyre Heros Bold Makie", halign = :center)
    if !isnothing(subtitle)
        Label(titlelayout[2, 1], subtitle, fontsize = 10)
    end
    rowgap!(titlelayout, 0)
    # rowgap!(f.layout, 0)
    # colgap!(f.layout, -50)
end

function hideall!(ax)
    tightlimits!(ax)
    hidedecorations!(ax)
    hidespines!(ax)
end

function axlims(x)
    min = minimum(x)
    max = maximum(x)
    return min - .1abs(min), max + .1abs(max)
end

function diagonal!(ax, x, y)
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    low = min(xmin, ymin)
    high = max(xmax, ymax)
    Makie.lines!(ax, [low, high], [low, high],
          color = :darkred, linewidth = 2)
end

function smoother!(ax, x, y, col = "red", crange = nothing, span = .5)
    us = range(extrema(x)...; step = .1)
    vs = Loess.predict(loess(x, y; span = span), us)
    if !isnothing(crange)
        f = Makie.lines!(ax, us, vs, color = col, colorrange = crange, linewidth = 4)
    else
        f = Makie.lines!(ax, us, vs, color = col, linewidth = 4)
    end
    return f
end

function calc_crange(x)
    m = maximum(abs.(x))
    return -m, m
end
