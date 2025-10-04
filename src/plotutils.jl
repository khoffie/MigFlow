getageyear(r) = r.mdl.meta.age, r.mdl.meta.year

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

function grid_position(k::Int, ncols::Int=2)
    i = div(k - 1, ncols) + 1   # row index
    j = mod(k - 1, ncols) + 1   # column index
    return (i, j)
end

function prettytitle!(f, title, row = 0)
    titlelayout = GridLayout(f[row, 1:end], tellwidth = false)
    Label(titlelayout[1, 1], title, fontsize = 15,
          font = "TeX Gyre Heros Bold Makie", halign = :center)
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
    return min - .1min, max + .1max
end
